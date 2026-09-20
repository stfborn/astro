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
# Verfügbare HiPS-Surveys
#
# WICHTIG zu "jwst": Es gibt KEIN allgemeines JWST-HiPS, das (wie DSS2/
# 2MASS/WISE) den ganzen Himmel abdeckt. "CDS/P/JWST/EPO" ist eine
# Sammlung von nur ~6 einzelnen, exakt gerahmten Outreach-Mosaiken
# (Carina-Nebel-Cliffs, Cartwheel-Galaxie, Southern Ring Nebula,
# Stephans Quintet, SMACS0723 Deep Field, First-Images-Übersicht) -
# Himmelsabdeckung ca. 0,00001 (0,001 %). Für praktisch alle Objekte
# liefert dieser Survey ein leeres/schwarzes Bild (siehe Warnung in
# download_images()).
SURVEYS = {
    "dss2":    ("CDS/P/DSS2/color",              "DSS2"),
    "hubble":  ("CDS/P/HST/color",               "HST"),
    "jwst":    ("CDS/P/JWST/EPO",                 "JWST"),
    "2mass":   ("CDS/P/2MASS/color",              "2MASS"),
    "wise":    ("CDS/P/WISE/W3",                  "WISE-W3"),
}
DEFAULT_SURVEY = "dss2"
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
    "gemini":      ["Castor", "Pollux", "Alhena", "Mebsuta", "Mekbuda", "Propus", "Tejat", "Wasat", "Alzirr"],
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
#
# WICHTIG: Simbad-Votable-Felder haben unterschiedliche Kardinalität!
#
#   - "sicherheitsrelevante" Felder wie V, B, sp_type, plx_value, pmra,
#     pmdec, otype liefern höchstens 1 Zeile pro Objekt (1:1 bzw. 1:0).
#   - Felder wie mesdiameter, mesrot, mesfe_h, rvz_radvel sind 1:n-
#     Messtabellen (mehrere publizierte Messwerte pro Objekt möglich).
#     Werden sie mit query_objects() für MEHRERE Namen gleichzeitig
#     abgefragt, entsteht ein Cross-Join: aus 3 angefragten Namen können
#     plötzlich >1000 Zeilen werden – und die naive Zuordnung
#     "Zeile i == Name i" liefert dann für falsche Sterne falsche Daten!
#
# Deshalb zwei getrennte Simbad-Instanzen:
#   SIM_BULK   -> nur 1:1-Felder, sicher für query_objects() (Konstellationen)
#   SIM_DETAIL -> zusätzliche 1:n-Felder, NUR für Einzelobjekt-Abfragen
#                 (query_object()). Da mehrere solcher 1:n-Tabellen
#                 gleichzeitig verknüpft werden, entstehen pro Objekt oft
#                 hunderte Kombinationszeilen OHNE garantierte Sortierung.
#                 Es wird daher NICHT einfach Zeile 0 genommen, sondern für
#                 jede Spalte einzeln über alle Zeilen der erste gültige Wert
#                 gesucht (siehe _first_valid_per_column) - so gehen keine
#                 vorhandenen Messwerte durch eine ungünstige Zeilen-
#                 Kombination verloren. Eine Verwechslung mit einem ANDEREN
#                 Stern ist ausgeschlossen, da hier immer nur ein Name
#                 abgefragt wird.
#
# Zusätzlich: SIM_DETAIL enthält bewusst KEINE V/B-Felder. Diese lösen
# einen Inner-Join auf die Fluss-Tabelle aus – Objekte ohne katalogisierte
# Breitbandphotometrie (z.B. Emissionsnebel wie M42) würden dadurch
# komplett aus dem Ergebnis herausfallen ("keine Daten gefunden"), obwohl
# Simbad durchaus Basisdaten (Objekttyp, Koordinaten, ...) dafür hat.
#
def _make_simbad_bulk():
    sim = Simbad()
    sim.add_votable_fields("V", "B", "sp_type", "plx_value", "pmra", "pmdec", "otype")
    return sim

def _make_simbad_detail():
    sim = Simbad()
    sim.add_votable_fields(
        "sp_type", "plx_value", "pmra", "pmdec", "otype",
        "mesdiameter", "mesrot", "mesfe_h", "rvz_radvel",
    )
    return sim

SIM_BULK   = _make_simbad_bulk()
SIM_DETAIL = _make_simbad_detail()

# ============================================================ #
# SIMBAD-STATISTIK: Zeile(n) in ein einfaches Dict umwandeln
# ============================================================ #
def _row_to_dict(row, colnames):
    """Wandelt eine Astropy-Tabellenzeile in ein {spalte: wert}-Dict um.
    Maskierte/fehlende Werte werden zu None."""
    out = {}
    for col in colnames:
        v = row[col]
        try:
            if hasattr(v, "mask") and v.mask:
                out[col] = None
                continue
        except (TypeError, ValueError):
            pass
        out[col] = v
    return out

def _is_empty(v):
    """Prüft, ob ein Simbad-Zellwert als 'leer/keine Daten' zu werten ist."""
    if v is None:
        return True
    try:
        if hasattr(v, "mask") and v.mask:
            return True
    except (TypeError, ValueError):
        pass
    s = str(v).strip()
    return s in ("", "--", "nan", "None")

def _first_valid_per_column(result):
    """
    Sucht für JEDE Spalte einzeln über ALLE Zeilen von `result` den ersten
    gültigen (nicht-leeren) Wert.

    Hintergrund: Werden mehrere 1:n-Messtabellen (mesrot, mesfe_h,
    mesdiameter, rvz_radvel) gemeinsam abgefragt, entsteht pro Objekt ein
    Cross-Join mit potenziell hunderten Kombinationszeilen - ohne
    garantierte Sortierung. Ein simples result[0] würde daher u.U. eine
    Kombination erwischen, in der z.B. die Metallizität zufällig leer ist,
    obwohl eine andere (ebenso gültige) Zeile einen echten Messwert hätte.
    Durch spaltenweises Suchen geht kein vorhandener Messwert verloren.
    """
    colnames = result.colnames
    out = {col: None for col in colnames}
    remaining = set(colnames)
    for row in result:
        if not remaining:
            break
        found = set()
        for col in remaining:
            v = row[col]
            if not _is_empty(v):
                out[col] = v
                found.add(col)
        remaining -= found
    return out

# ============================================================ #
# SIMBAD-STATISTIK: Einzelobjekt ausgeben
# ============================================================ #
def _print_simbad_row(name, data):
    """Gibt die Simbad-Statistik für ein Objekt aus.
    `data` ist ein {spalte: wert}-Dict, ggf. aus mehreren Abfragen
    (Basis + Detailfelder) zusammengeführt."""

    def val(col):
        v = data.get(col)
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
# SIMBAD: Mehrere Objekte auf einmal abfragen (nur sichere 1:1-Felder!)
# ============================================================ #
def query_simbad_bulk(names):
    """
    Fragt Simbad für eine Liste von Namen in einer einzigen Anfrage ab.
    Gibt ein Dict {name: {spalte: wert}} zurück.

    Verwendet ausschließlich SIM_BULK (Felder mit Kardinalität <= 1), damit
    result[i] garantiert zum i-ten Namen passt. Sollte die Zeilenzahl trotzdem
    nicht zur Namenszahl passen (z.B. durch zukünftige Feldänderungen), wird
    das Ergebnis sicherheitshalber verworfen statt falsch zugeordnet.
    """
    try:
        result = SIM_BULK.query_objects(names)
    except Exception as e:
        tprint(f"  [Simbad] Bulk-Abfrage fehlgeschlagen: {e}")
        return {}

    if result is None or len(result) == 0:
        return {}

    if len(result) != len(names):
        tprint(f"  [Simbad] Warnung: Bulk-Abfrage lieferte {len(result)} Zeilen "
               f"für {len(names)} Namen (Cross-Join?) - Ergebnis wird verworfen, "
               f"um Falschzuordnungen zu vermeiden.")
        return {}

    out = {}
    colnames = result.colnames
    for i, name in enumerate(names):
        out[name] = _row_to_dict(result[i], colnames)
    return out

# ============================================================ #
# SIMBAD: Detailfelder (1:n-Messtabellen) für EIN Objekt nachladen
# ============================================================ #
def query_simbad_detail(name):
    """
    Holt zusätzliche 1:n-Messfelder (Durchmesser, Rotation, Metallizität,
    Radialgeschwindigkeit) für ein einzelnes Objekt. Da diese Felder pro
    Objekt viele Kombinationszeilen liefern können, wird über
    _first_valid_per_column() für jede Spalte einzeln der erste gültige
    Wert über alle Zeilen gesucht - so gehen keine vorhandenen Messwerte
    verloren. Eine Verwechslung mit einem anderen Objekt ist ausgeschlossen,
    da hier immer genau ein Name abgefragt wird.
    """
    try:
        result = SIM_DETAIL.query_object(name)
    except Exception:
        return {}
    if result is None or len(result) == 0:
        return {}
    return _first_valid_per_column(result)

# ============================================================ #
# SIMBAD: Einzelobjekt (Basis- + Detailfelder kombiniert)
# ============================================================ #
def query_simbad_single(name):
    data = {}
    try:
        result = SIM_BULK.query_object(name)
        if result is not None and len(result) > 0:
            data.update(_row_to_dict(result[0], result.colnames))
    except Exception as e:
        tprint(f"  [Simbad] Abfrage fehlgeschlagen: {e}")

    data.update(query_simbad_detail(name))

    if not data:
        tprint("  [Simbad] Keine Daten gefunden.")
        return
    _print_simbad_row(name, data)

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
def _check_png_coverage(content):
    """
    Prüft, ob ein heruntergeladenes PNG praktisch leer/schwarz ist (keine
    tatsächliche Bilddaten an der angefragten Position). Kommt vor allem
    bei Surveys mit lückenhafter Abdeckung vor (Hubble ~0,18 %, JWST
    ~0,001 % des Himmels) - der HiPS-Dienst liefert dann trotzdem ein
    technisch gültiges, aber inhaltsloses Bild statt eines Fehlers.

    Gibt True zurück, wenn das Bild als "leer" zu werten ist, sonst False.
    Bei Fehlern beim Öffnen (z.B. Pillow fehlt) wird konservativ False
    zurückgegeben - dann eben ohne Warnung.
    """
    try:
        from PIL import Image
        import numpy as np
        import io
        img = Image.open(io.BytesIO(content)).convert("L")
        arr = np.asarray(img, dtype=np.float64)
        # Mittelwert statt Maximalwert, da Kachel-/Rand-Antialiasing auch bei
        # inhaltsleeren Bildern vereinzelt helle Randpixel erzeugen kann.
        # Empirisch: echte Aufnahmen liegen deutlich über mean=10 (z.B. DSS2
        # M1 ≈ 8.9, Hubble M104-Ausschnitt ≈ 70, Vega-Punktquelle ≈ 119),
        # leere HiPS-Kacheln (keine Abdeckung) liegen bei mean < 1.5.
        return arr.mean() < 2.0
    except Exception:
        return False

def download_images(coord, label, fov=DEFAULT_FOV, formats=("png", "fits"),
                     surveys=("dss2",), projection="SIN"):
    """
    Lädt Bilder für die angegebenen Surveys/Formate herunter.

    WICHTIG: Ganzhimmels-Projektionen (AIT = Aitoff, MOL = Mollweide) zeigen
    per Definition immer die komplette Himmelskugel. Das übergebene `fov`
    wird für diese Projektionen daher IGNORIERT und stattdessen automatisch
    auf 360° gesetzt - unabhängig davon, was per --fov angegeben wurde.
    """
    safe_name = label.replace(" ", "_").replace("/", "_")

    proj = projection.upper()
    effective_fov = 360 if proj in ("AIT", "MOL") else fov
    suffix = "" if proj == "SIN" else f"_{proj}"

    base_params = {
        "width": SIZE, "height": SIZE,
        "ra": coord.ra.deg, "dec": coord.dec.deg,
        "coordsys": "icrs", "fov": effective_fov,
        "projection": proj, "min_cut": "0.5%", "max_cut": "99.5%", "stretch": "linear",
    }
    fov_hinweis = f"{effective_fov}°" if effective_fov == fov else f"{effective_fov}° (Ganzhimmel, --fov ignoriert)"
    tprint(f"\nBild-Download: Sichtfeld = {fov_hinweis}, Projektion: {proj}, "
           f"Surveys: {', '.join(surveys).upper()}...")
    for survey_key in surveys:
        hips_id, survey_label = SURVEYS[survey_key]
        for fmt in formats:
            filename = OUTPUT_DIR / f"{safe_name}_{survey_label}{suffix}.{fmt}"
            params   = {**base_params, "hips": hips_id, "format": fmt}
            for url in (HIPS_URL, HIPS_URL_MIRROR):
                try:
                    # Timeout als (Connect, Read)-Tupel: Ist ein Server
                    # komplett down (TCP-Handshake antwortet nicht), wird das
                    # schnell erkannt (8s) statt erst nach 120s auf den
                    # Mirror auszuweichen. Für langsame, aber laufende
                    # Downloads bleibt weiterhin viel Zeit (100s Read).
                    r = requests.get(url, params=params, timeout=(8, 100))
                    r.raise_for_status()
                    filename.write_bytes(r.content)
                    hinweis = ""
                    if fmt == "png" and _check_png_coverage(r.content):
                        hinweis = ("  ⚠ Bild wirkt leer/schwarz - vermutlich keine "
                                   f"{survey_label}-Abdeckung an dieser Position")
                    tprint(f"  -> [{survey_label}] {filename} ({len(r.content):,} Bytes){hinweis}")
                    break
                except Exception as e:
                    tprint(f"  -> [{survey_label}] Download über {url} fehlgeschlagen: {e}")

# ============================================================ #
# EINZELOBJEKT verarbeiten
# ============================================================ #
def process_single_object(name, fov, formats=("png",), surveys=("dss2",), projection="SIN"):
    tprint("=" * 60)
    tprint(f" Objekt: {name}  |  Sichtfeld: {fov}°")
    tprint("=" * 60)

    coord, source = resolve_object(name)
    if coord is None:
        tprint(f"  -> FEHLER: '{name}' konnte nicht identifiziert werden.\n")
        return False

    tprint(f"  -> Aufgelöst via {source}: RA {coord.ra.deg:.6f}°  Dec {coord.dec.deg:.6f}°")
    query_simbad_single(name)
    download_images(coord, name, fov=fov, formats=formats, surveys=surveys, projection=projection)
    tprint()
    return True

# ============================================================ #
# KONSTELLATION: parallel verarbeiten
# ============================================================ #
def resolve_and_report_star(name, simbad_cache):
    """
    Löst einen Stern auf und gibt (falls vorhanden) seine Simbad-Statistik aus.
    Wird sowohl von process_constellation() als auch von constellation_gesamt()
    verwendet, damit die Detailausgabe in beiden Modi identisch ist.

    Basisdaten kommen aus dem Bulk-Cache (sicher zugeordnet), Detailfelder
    (Durchmesser/Rotation/Metallizität/RV) werden pro Stern einzeln
    nachgeladen - das läuft ohnehin bereits in einem eigenen Thread und
    kann daher nicht mit einem anderen Stern verwechselt werden.

    Gibt das aufgelöste SkyCoord zurück (oder None bei Fehlschlag).
    """
    coord, source = resolve_object(name)
    if coord is None:
        tprint(f"  [!] '{name}' konnte nicht identifiziert werden.")
        return None

    tprint(f"  -> {name}: aufgelöst via {source} "
           f"(RA {coord.ra.deg:.4f}°  Dec {coord.dec.deg:.4f}°)")

    data = dict(simbad_cache.get(name, {}))
    data.update(query_simbad_detail(name))
    if data:
        _print_simbad_row(name, data)

    return coord

def process_star_parallel(name, fov, formats, surveys, simbad_cache, projection="SIN"):
    coord = resolve_and_report_star(name, simbad_cache)
    if coord is None:
        return False

    download_images(coord, name, fov=fov, formats=formats, surveys=surveys, projection=projection)
    return True

def process_constellation(stars, fov, formats, surveys, projection="SIN"):
    tprint(f"\n  Frage Simbad für {len(stars)} Sterne auf einmal ab...")
    simbad_cache = query_simbad_bulk(stars)
    tprint(f"  -> {len(simbad_cache)}/{len(stars)} Einträge erhalten.\n")

    successes = 0
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {
            executor.submit(process_star_parallel, star, fov, formats, surveys, simbad_cache, projection): star
            for star in stars
        }
        for future in as_completed(futures):
            if future.result():
                successes += 1
    return successes

def constellation_gesamt(constellation_name, stars, fov, formats, surveys, projection="SIN"):
    """
    Löst alle Hauptsterne auf (inkl. Simbad-Detailausgabe je Stern), berechnet
    den Mittelpunkt der Konstellation und lädt ein einzelnes Übersichtsbild
    herunter.
    """
    tprint(f"\n  Frage Simbad für {len(stars)} Sterne auf einmal ab...")
    simbad_cache = query_simbad_bulk(stars)
    tprint(f"  -> {len(simbad_cache)}/{len(stars)} Einträge erhalten.\n")

    tprint(f"  Löse Koordinaten für {len(stars)} Hauptsterne auf...")

    coords = []
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {
            executor.submit(resolve_and_report_star, star, simbad_cache): star
            for star in stars
        }
        for future in as_completed(futures):
            coord = future.result()
            if coord:
                coords.append(coord)

    if not coords:
        tprint("FEHLER: Keine Koordinaten aufgelöst.")
        return False

    # Mittelpunkt berechnen (kartesisch mitteln, dann zurück in RA/Dec)
    all_coords = SkyCoord(coords)
    center = all_coords.cartesian.mean()
    center_coord = SkyCoord(center, representation_type="cartesian").represent_as("unitspherical")
    center_skycoord = SkyCoord(ra=center_coord.lon, dec=center_coord.lat, frame="icrs")

    tprint(f"\n  Mittelpunkt: RA {center_skycoord.ra.deg:.4f}°  Dec {center_skycoord.dec.deg:.4f}°")
    tprint(f"  Lade Gesamtbild (FOV = {fov}°)...")

    label = f"{constellation_name.title()}_gesamt"
    download_images(center_skycoord, label, fov=fov, formats=formats, surveys=surveys, projection=projection)
    return True


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
        "--constellation", metavar="NAME", nargs="+",
        help="Lädt alle Hauptsterne einer Konstellation, z.B. Orion oder "
             "Ursa Major (mehrere Wörter auch ohne Anführungszeichen möglich). "
             f"Verfügbar: {', '.join(sorted(set(CONSTELLATION_MAIN_STARS)))}",
    )
    parser.add_argument(
        "--format", choices=["png", "fits", "both"], default=None,
        help="Ausgabeformat: png, fits oder both. "
             "Default: both beim Einzelobjekt, png bei --constellation.",
    )
    parser.add_argument(
        "--survey",
        nargs="+",
        choices=list(SURVEYS.keys()) + ["all"],
        default=[DEFAULT_SURVEY],
        metavar="NAME",
        help=f"HiPS-Survey(s): {', '.join(SURVEYS.keys())}, all. "
             f"Mehrere möglich: --survey dss2 hubble. Default: {DEFAULT_SURVEY}",
    )
    parser.add_argument(
        "--gesamt",
        action="store_true",
        help="Nur bei --constellation: lädt ein einziges Übersichtsbild "
             "zentriert auf den Mittelpunkt der Konstellation.",
    )
    parser.add_argument(
        "--projection",
        choices=["SIN", "TAN", "ARC", "ZEA", "STG", "CAR", "AIT", "MOL"],
        default="SIN", type=str.upper,
        help="WCS-Projektion für das Bild (Default: SIN). AIT (Aitoff) und "
             "MOL (Mollweide) sind Ganzhimmels-Projektionen: bei diesen wird "
             "automatisch das komplette 360°-Sichtfeld verwendet, "
             "unabhängig vom --fov-Wert.",
    )
    args = parser.parse_args()

    # Survey-Liste auflösen
    if "all" in args.survey:
        surveys = list(SURVEYS.keys())
    else:
        surveys = args.survey

    def resolve_formats(default):
        if args.format == "png":  return ("png",)
        if args.format == "fits": return ("fits",)
        if args.format == "both": return ("png", "fits")
        return default

    if args.constellation:
        constellation_name = " ".join(args.constellation)
        key   = constellation_name.strip().lower()
        stars = CONSTELLATION_MAIN_STARS.get(key)
        if stars is None:
            print(f"FEHLER: '{constellation_name}' nicht in der Liste.")
            print(f"Verfügbar: {', '.join(sorted(set(CONSTELLATION_MAIN_STARS)))}")
            sys.exit(1)

        formats = resolve_formats(("png",))
        print("=" * 60)
        print(f" Konstellation: {constellation_name.title()}  "
              f"({len(stars)} Hauptsterne, Sichtfeld = {args.fov}°)")
        print("=" * 60)

        if args.gesamt:
            constellation_gesamt(constellation_name, stars, args.fov, formats, surveys, projection=args.projection)
        else:
            print(f" {MAX_WORKERS} parallele Threads")
            successes = process_constellation(stars, args.fov, formats, surveys, projection=args.projection)
            print("=" * 60)
            print(f" Fertig. {successes}/{len(stars)} Sterne erfolgreich heruntergeladen.")

        print("=" * 60)
        return

    if not args.object:
        parser.print_help()
        sys.exit(1)

    search_query = " ".join(args.object)
    process_single_object(search_query, args.fov, formats=resolve_formats(("png", "fits")),
                           surveys=surveys, projection=args.projection)

if __name__ == "__main__":
    main()
