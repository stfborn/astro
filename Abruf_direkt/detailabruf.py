#!/usr/bin/env python3

import sys
import argparse
import math
import threading
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad

# ============================================================ #
# EINSTELLUNGEN
# ============================================================ #
DEFAULT_FOV  = 3.0   # Sichtfeld in Grad (Default, per --fov überschreibbar)
SIZE         = 2048  # Bildgröße in Pixeln
MAX_WORKERS  = 5     # Parallele Threads für Konstellations-Downloads

# Offizielle, stabile CDS-Dienste
HIPS_URL        = "https://alasky.cds.unistra.fr/hips-image-services/hips2fits"
HIPS_URL_MIRROR = "https://alaskybis.cds.unistra.fr/hips-image-services/hips2fits"
HIPS_SURVEY     = "CDS/P/DSS2/color"
OUTPUT_DIR      = Path(".")

# Spaltennamen, unter denen VizieR-Kataloge Koordinaten üblicherweise führen
RA_COLUMNS  = ("RAJ2000", "RA_ICRS", "_RAJ2000", "RAICRS")
DEC_COLUMNS = ("DEJ2000", "DE_ICRS", "_DEJ2000", "DEICRS")

# Lock für saubere Konsolenausgabe bei parallelen Threads
PRINT_LOCK = threading.Lock()

def tprint(*args, **kwargs):
    """Thread-sicheres print()."""
    with PRINT_LOCK:
        print(*args, **kwargs)

# Hauptsterne (Bayer-Bezeichnung) gängiger Konstellationen
CONSTELLATION_MAIN_STARS = {
    "orion":       ["Betelgeuse", "Rigel", "Bellatrix", "Mintaka", "Alnilam", "Alnitak", "Saiph"],
    "ursa major":  ["Dubhe", "Merak", "Phecda", "Megrez", "Alioth", "Mizar", "Alkaid"],
    "cassiopeia":  ["Schedar", "Caph", "Gamma Cassiopeiae", "Ruchbah", "Segin"],
    "cygnus":      ["Deneb", "Albireo", "Sadr", "Gienah", "Delta Cygni"],
    "scorpius":    ["Antares", "Shaula", "Sargas", "Dschubba", "Graffias"],
    "leo":         ["Regulus", "Denebola", "Algieba", "Zosma", "Chertan"],
    "taurus":      ["Aldebaran", "Elnath", "Alcyone", "Zeta Tauri"],
    "gemini":      ["Castor", "Pollux", "Alhena", "Mebsuta"],
    "canis major": ["Sirius", "Mirzam", "Wezen", "Adhara"],
    "ursa minor":  ["Polaris", "Kochab", "Pherkad"],
    "lyra":        ["Vega", "Sheliak", "Sulafat"],
    "aquila":      ["Altair", "Tarazed", "Alshain"],
    "boötes":      ["Arcturus", "Izar", "Muphrid", "Seginus"],
    "bootes":      ["Arcturus", "Izar", "Muphrid", "Seginus"],
    "andromeda":   ["Alpheratz", "Mirach", "Almach"],
    "pegasus":     ["Enif", "Scheat", "Markab", "Algenib"],
    "perseus":     ["Mirfak", "Algol", "Atik"],
    "auriga":      ["Capella", "Menkalinan", "Almaaz"],
    "canis minor": ["Procyon", "Gomeisa"],
    "sagittarius": ["Kaus Australis", "Nunki", "Ascella"],
}

# ============================================================ #
# SIMBAD-KONFIGURATION (einmalig erstellen)
# ============================================================ #
def _make_simbad():
    sim = Simbad()
    sim.add_votable_fields(
        "V", "B", "sp_type", "plx_value", "pmra", "pmdec",
        "otype", "mesdiameter", "mesrot", "mesfe_h", "rvz_radvel",
    )
    return sim

SIM = _make_simbad()

# ============================================================ #
# SIMBAD-STATISTIK: Einzelobjekt ausgeben
# ============================================================ #
def _print_simbad_row(name, row, colnames):
    """Gibt die Simbad-Statistik für eine Zeile aus."""

    def val(col):
        if col not in colnames:
            return None
        v = row[col]
        try:
            if hasattr(v, 'mask') and v.mask:
                return None
        except (TypeError, ValueError):
            pass
        if v is None:
            return None
        s = str(v).strip()
        return s if s not in ("", "--", "nan", "None") else None

    lines = [f"\n{'─' * 60}", f" SIMBAD-STATISTIK: {name}", "─" * 60]

    main_id = val("main_id")
    otype   = val("otype")
    if main_id: lines.append(f"  Hauptbezeichnung : {main_id}")
    if otype:   lines.append(f"  Objekttyp        : {otype}")

    ra  = val("ra")
    dec = val("dec")
    if ra and dec:
        lines.append(f"  RA / Dec         : {ra}  /  {dec}")

    sp = val("sp_type")
    if sp:
        lines.append(f"  Spektraltyp      : {sp}")

    vmag = val("V")
    bmag = val("B")
    if vmag: lines.append(f"  Helligkeit (V)   : {float(vmag):.2f} mag")
    if bmag: lines.append(f"  Helligkeit (B)   : {float(bmag):.2f} mag")

    plx = val("plx_value")
    if plx:
        plx_f = float(plx)
        if plx_f > 0:
            dist_pc = 1000.0 / plx_f
            dist_ly = dist_pc * 3.26156
            lines.append(f"  Parallaxe        : {plx_f:.4f} mas")
            lines.append(f"  Entfernung       : {dist_ly:,.1f} Lichtjahre  ({dist_pc:,.1f} Parsec)")
            if vmag:
                abs_mag = float(vmag) - 5 * math.log10(dist_pc / 10)
                lines.append(f"  Abs. Helligkeit  : {abs_mag:.2f} mag")

    pm_ra  = val("pmra")
    pm_dec = val("pmdec")
    if pm_ra and pm_dec:
        pm_total = math.hypot(float(pm_ra), float(pm_dec))
        lines.append(f"  Eigenbewegung    : {float(pm_ra):.2f} mas/yr (RA)  "
                     f"{float(pm_dec):.2f} mas/yr (Dec)  → gesamt {pm_total:.2f} mas/yr")

    rv = val("rvz_radvel")
    if rv:
        rv_f = float(rv)
        richtung = "auf Baryzentrum zu" if rv_f < 0 else "vom Baryzentrum weg"
        lines.append(f"  Radialgeschw.    : {rv_f:+.1f} km/s ({richtung}, barizentrisch)")

    vsini = val("mesrot.vsini")
    if vsini:
        lines.append(f"  v·sin(i)         : {float(vsini):.1f} km/s")

    feh = val("mesfe_h.fe_h")
    if feh:
        lines.append(f"  Metallizität     : [Fe/H] = {float(feh):.2f}")

    teff = val("mesfe_h.teff")
    if teff:
        lines.append(f"  Effektivtemp.    : {int(float(teff))} K")

    diam      = val("mesdiameter.diameter")
    diam_unit = val("mesdiameter.unit")
    if diam:
        lines.append(f"  Winkelausdehnung : {float(diam):.2f} {diam_unit or 'arcsec'}")

    lines.append("─" * 60)
    tprint("\n".join(lines))

# ============================================================ #
# SIMBAD: Mehrere Objekte auf einmal abfragen
# ============================================================ #
def query_simbad_bulk(names):
    """
    Fragt Simbad für eine Liste von Namen in einer einzigen Anfrage ab.
    Gibt ein Dict {name: (row, colnames)} zurück.
    """
    try:
        result = SIM.query_objects(names)
    except Exception as e:
        tprint(f"  [Simbad] Bulk-Abfrage fehlgeschlagen: {e}")
        return {}

    if result is None or len(result) == 0:
        return {}

    out = {}
    colnames = result.colnames
    for i, name in enumerate(names):
        if i < len(result):
            out[name] = (result[i], colnames)
    return out

# ============================================================ #
# SIMBAD: Einzelobjekt (für Einzelaufruf)
# ============================================================ #
def query_simbad_single(name):
    try:
        result = SIM.query_object(name)
    except Exception as e:
        tprint(f"  [Simbad] Abfrage fehlgeschlagen: {e}")
        return
    if result is None or len(result) == 0:
        tprint("  [Simbad] Keine Daten gefunden.")
        return
    _print_simbad_row(name, result[0], result.colnames)

# ============================================================ #
# OBJEKT-AUFLÖSUNG
# ============================================================ #
def resolve_via_sesame(name):
    try:
        coord = SkyCoord.from_name(name)
        return coord
    except Exception:
        return None

def resolve_via_vizier(name):
    try:
        result = Vizier(columns=["*"]).query_object(name)
        if not result:
            return None
        for table in result:
            ra_col  = next((c for c in RA_COLUMNS  if c in table.colnames), None)
            dec_col = next((c for c in DEC_COLUMNS if c in table.colnames), None)
            if ra_col and dec_col:
                return SkyCoord(
                    ra=float(table[ra_col][0]) * u.deg,
                    dec=float(table[dec_col][0]) * u.deg,
                    frame="icrs"
                )
    except Exception:
        pass
    return None

def resolve_object(name):
    coord = resolve_via_sesame(name)
    source = "Sesame"
    if coord is None:
        coord = resolve_via_vizier(name)
        source = "VizieR"
    return coord, source

# ============================================================ #
# BILD-DOWNLOAD
# ============================================================ #
def download_dss2_images(coord, label, fov=DEFAULT_FOV, formats=("png", "fits")):
    safe_name   = label.replace(" ", "_").replace("/", "_")
    base_params = {
        "hips": HIPS_SURVEY, "width": SIZE, "height": SIZE,
        "ra": coord.ra.deg, "dec": coord.dec.deg,
        "coordsys": "icrs", "fov": fov,
        "projection": "SIN", "min_cut": "0.5%", "max_cut": "99.5%", "stretch": "linear",
    }
    for fmt in formats:
        filename = OUTPUT_DIR / f"{safe_name}_DSS2.{fmt}"
        params   = {**base_params, "format": fmt}
        for url in (HIPS_URL, HIPS_URL_MIRROR):
            try:
                r = requests.get(url, params=params, timeout=120)
                r.raise_for_status()
                filename.write_bytes(r.content)
                tprint(f"  -> Gespeichert: {filename} ({len(r.content):,} Bytes)")
                break
            except Exception as e:
                tprint(f"  -> Download über {url} fehlgeschlagen: {e}")

# ============================================================ #
# EINZELOBJEKT verarbeiten
# ============================================================ #
def process_single_object(name, fov, formats=("png",)):
    tprint("=" * 60)
    tprint(f" Objekt: {name}  |  Sichtfeld: {fov}°")
    tprint("=" * 60)

    coord, source = resolve_object(name)
    if coord is None:
        tprint(f"  -> FEHLER: '{name}' konnte nicht identifiziert werden.\n")
        return False

    tprint(f"  -> Aufgelöst via {source}: RA {coord.ra.deg:.6f}°  Dec {coord.dec.deg:.6f}°")
    query_simbad_single(name)
    tprint(f"\nDSS2 API: Starte Bild-Download (Sichtfeld = {fov} Grad)...")
    download_dss2_images(coord, name, fov=fov, formats=formats)
    tprint()
    return True

# ============================================================ #
# KONSTELLATION: parallel verarbeiten
# ============================================================ #
def process_star_parallel(name, fov, formats, simbad_cache):
    """Wird im Thread ausgeführt: Koordinaten auflösen + Bild laden + Stats ausgeben."""
    coord, source = resolve_object(name)
    if coord is None:
        tprint(f"  [!] '{name}' konnte nicht identifiziert werden.")
        return False

    tprint(f"  -> {name}: aufgelöst via {source} "
           f"(RA {coord.ra.deg:.4f}°  Dec {coord.dec.deg:.4f}°)")

    # Simbad-Statistik aus dem vorbereiteten Cache ausgeben
    if name in simbad_cache:
        row, colnames = simbad_cache[name]
        _print_simbad_row(name, row, colnames)

    download_dss2_images(coord, name, fov=fov, formats=formats)
    return True

def process_constellation(stars, fov, formats):
    tprint(f"\n  Frage Simbad für {len(stars)} Sterne auf einmal ab...")
    simbad_cache = query_simbad_bulk(stars)
    tprint(f"  -> {len(simbad_cache)}/{len(stars)} Einträge erhalten.\n")

    successes = 0
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {
            executor.submit(process_star_parallel, star, fov, formats, simbad_cache): star
            for star in stars
        }
        for future in as_completed(futures):
            if future.result():
                successes += 1
    return successes

# ============================================================ #
# HAUPTPROGRAMM
# ============================================================ #
def main():
    parser = argparse.ArgumentParser(
        description="Lädt DSS2-Bilder für Himmelsobjekte herunter und zeigt Simbad-Statistiken."
    )
    parser.add_argument(
        "object", nargs="*",
        help="Objektname, z.B. M42, NGC7000, Vega, HD39801 (entfällt bei --constellation)"
    )
    parser.add_argument(
        "--fov", type=float, default=DEFAULT_FOV,
        help=f"Sichtfeld in Grad (Default: {DEFAULT_FOV}). Für Sterne z.B. 0.1–0.5.",
    )
    parser.add_argument(
        "--constellation", metavar="NAME",
        help="Lädt alle Hauptsterne einer Konstellation, z.B. Orion. "
             f"Verfügbar: {', '.join(sorted(set(CONSTELLATION_MAIN_STARS)))}",
    )
    parser.add_argument(
        "--format", choices=["png", "fits", "both"], default=None,
        help="Ausgabeformat: png, fits oder both. "
             "Default: both beim Einzelobjekt, png bei --constellation.",
    )
    args = parser.parse_args()

    def resolve_formats(default):
        if args.format == "png":  return ("png",)
        if args.format == "fits": return ("fits",)
        if args.format == "both": return ("png", "fits")
        return default

    if args.constellation:
        key   = args.constellation.strip().lower()
        stars = CONSTELLATION_MAIN_STARS.get(key)
        if stars is None:
            print(f"FEHLER: '{args.constellation}' nicht in der Liste.")
            print(f"Verfügbar: {', '.join(sorted(set(CONSTELLATION_MAIN_STARS)))}")
            sys.exit(1)

        formats = resolve_formats(("png",))
        print("=" * 60)
        print(f" Konstellation: {args.constellation.title()}  "
              f"({len(stars)} Hauptsterne, Sichtfeld = {args.fov}°, "
              f"{MAX_WORKERS} parallele Threads)")
        print("=" * 60)

        successes = process_constellation(stars, args.fov, formats)

        print("=" * 60)
        print(f" Fertig. {successes}/{len(stars)} Sterne erfolgreich heruntergeladen.")
        print("=" * 60)
        return

    if not args.object:
        parser.print_help()
        sys.exit(1)

    search_query = " ".join(args.object)
    process_single_object(search_query, args.fov, formats=resolve_formats(("png", "fits")))

if __name__ == "__main__":
    main()
