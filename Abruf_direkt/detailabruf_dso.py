#!/usr/bin/env python3
"""
detailabruf_dso.py – Deep-Sky-Objekte (Galaxien, Nebel, Sternhaufen)

Variante von detailabruf_aktualisiert.py, exklusiv auf Deep-Sky-Objekte
zugeschnitten. Unterschiede zur Stern-Version:

  - SIMBAD-Felder auf DSO-Relevanz umgestellt: Winkelausdehnung
    (galdim_majaxis/minaxis/angle), Morphologie (morph_type),
    Radialgeschwindigkeit/Rotverschiebung (rvz_radvel/rvz_redshift)
    statt Parallaxe/Spektraltyp/Eigenbewegung.
  - V/B-Helligkeit wird separat abgefragt (löst bei DSOs ohne
    Breitbandphotometrie sonst einen Inner-Join aus, der das Objekt
    komplett aus dem Ergebnis verschwinden lässt - z.B. bei M42/NGC7000).
  - Konstellationsliste enthält bekannte Messier/NGC-Objekte statt
    Hauptsterne.
  - Default-Sichtfeld 1.0° (typische DSO-Größenordnung) statt 3.0°.
"""

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
DEFAULT_FOV  = 1.0   # Sichtfeld in Grad (Default, per --fov überschreibbar) - typische DSO-Größe
SIZE         = 2048  # Bildgröße in Pixeln
MAX_WORKERS  = 5     # Parallele Threads für Konstellations-Downloads

# Lichtgeschwindigkeit (für Hubble-Distanz aus Rotverschiebung, näherungsweise)
C_KM_S = 299792.458
H0_KM_S_MPC = 70.0  # Hubble-Konstante (Näherungswert)

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
# Himmelsabdeckung ca. 0,00001 (0,001 %). Selbst für die o.g. Zielobjekte
# (z.B. NGC 3372/Carina, Stephan's Quintet) trifft eine Abfrage über die
# SIMBAD-Katalogkoordinaten des Gesamtobjekts die winzige, exakt gerahmte
# JWST-Kachel NICHT zuverlässig, da deren Zentrum von der SIMBAD-Position
# abweicht (z.B. NIRCam-Ausschnitt "Cosmic Cliffs" liegt in NGC 3324, einem
# kleinen Teilbereich von NGC 3372, nicht am SIMBAD-Zentrum von NGC 3372).
# Für praktisch jede Abfrage liefert dieser Survey daher ein leeres/
# schwarzes Bild (siehe Warnung in download_images()) - "jwst" eignet sich
# aktuell nicht für allgemeine Objektabfragen.
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

# Bekannte Deep-Sky-Objekte (Messier/NGC) gängiger Konstellationen
CONSTELLATION_MAIN_STARS = {
    "orion":       ["M42", "M43", "M78", "NGC 2024", "IC 434", "M1"],
    "ursa major":  ["M81", "M82", "M97", "M101", "M108", "M109"],
    "cassiopeia":  ["M52", "M103", "NGC 457", "NGC 663"],
    "cygnus":      ["NGC 7000", "M39", "NGC 6960", "NGC 6992", "NGC 6826"],
    "scorpius":    ["M4", "M6", "M7", "M80", "NGC 6231"],
    "leo":         ["M65", "M66", "M95", "M96", "M105"],
    "taurus":      ["M1", "M45"],
    "gemini":      ["M35", "NGC 2392"],
    "canis major": ["M41", "NGC 2359"],
    "ursa minor":  ["NGC 6217"],
    "lyra":        ["M57", "M56"],
    "aquila":      ["NGC 6709", "NGC 6741"],
    "boötes":      ["NGC 5466"],
    "bootes":      ["NGC 5466"],
    "andromeda":   ["M31", "M32", "M110"],
    "pegasus":     ["M15", "NGC 7331", "Stephan's Quintet"],
    "perseus":     ["M34", "M76", "NGC 869", "NGC 884"],
    "auriga":      ["M36", "M37", "M38"],
    "canis minor": [],
    "sagittarius": ["M8", "M20", "M17", "M22", "M23", "M28"],
    "vulpecula":   ["M27"],
    "hercules":    ["M13", "M92"],
    "monoceros":   ["NGC 2244", "M50"],
    "sagitta":     ["M71"],
    "vela":        ["NGC 3132"],
    "carina":      ["NGC 3372", "NGC 3532"],
    "puppis":      ["M46", "M47", "M93"],
    "canes venatici": ["M51", "M63", "M94", "M3"],
    "coma berenices": ["M64", "M85", "M88", "M100", "Coma Cluster"],
    "virgo":       ["M49", "M58", "M59", "M60", "M84", "M86", "M87", "M104"],
    "sculptor":    ["NGC 253", "NGC 55"],
    "triangulum":  ["M33"],
    "aquarius":    ["M2", "M72", "M73", "NGC 7009"],
    "capricornus": ["M30"],
    "ophiuchus":   ["M9", "M10", "M12", "M14", "M19", "M62", "M107"],
    "serpens":     ["M5", "M16"],
    "cancer":      ["M44", "M67"],
    "fornax":      ["NGC 1365", "Fornax Cluster"],
}

# Bezeichnung "CONSTELLATION_DSOS" als klarerer Alias (gleiche Daten)
CONSTELLATION_DSOS = CONSTELLATION_MAIN_STARS

# ============================================================ #
# SIMBAD-KONFIGURATION (einmalig erstellen) - DSO-VARIANTE
# ============================================================ #
#
# WICHTIG: Simbad-Votable-Felder haben unterschiedliche Kardinalität!
#
#   - otype, morphtype (morph_type/-qual/-bibcode), dimensions
#     (galdim_majaxis/minaxis/angle/qual/bibcode) und velocity
#     (rvz_radvel/rvz_redshift/...) sind alle "column of basic",
#     liefern also höchstens 1 Zeile pro Objekt (1:1 bzw. 1:0) und
#     sind damit SICHER für query_objects() (Bulk-Abfrage mehrerer
#     Namen gleichzeitig).
#   - mesDistance ist dagegen eine 1:n-Messtabelle (mehrere publizierte
#     Entfernungsbestimmungen pro Objekt, verschiedene Methoden/Einheiten).
#     Wird sie mit query_objects() für MEHRERE Namen gleichzeitig
#     abgefragt, entsteht ein Cross-Join - die naive Zuordnung
#     "Zeile i == Name i" würde dann für falsche Objekte falsche Daten
#     liefern! Deshalb NUR per query_object() für ein einzelnes Objekt.
#
# Deshalb zwei getrennte Simbad-Instanzen:
#   SIM_BULK   -> nur 1:1-Felder, sicher für query_objects() (Konstellationen)
#   SIM_DETAIL -> zusätzliche 1:n-Felder (mesDistance), NUR für
#                 Einzelobjekt-Abfragen (query_object()).
#
# Zusätzlich: WEDER SIM_BULK NOCH SIM_DETAIL enthalten V/B-Felder!
# Diese lösen einen Inner-Join auf die Fluss-Tabelle aus - Objekte ohne
# katalogisierte Breitbandphotometrie (z.B. Emissionsnebel wie M42,
# NGC 7000) würden dadurch komplett aus dem Ergebnis herausfallen
# ("keine Daten gefunden"), obwohl Simbad durchaus Basisdaten
# (Objekttyp, Größe, Koordinaten, ...) dafür hat. V/B wird daher über
# eine dritte, isolierte Simbad-Instanz (SIM_PHOTOMETRY) abgefragt und
# nur bei Erfolg zusammengeführt - fehlt sie, bleiben alle anderen
# Felder trotzdem erhalten.
#
def _make_simbad_bulk():
    sim = Simbad()
    sim.add_votable_fields("otype", "morphtype", "dimensions", "velocity")
    return sim

def _make_simbad_detail():
    sim = Simbad()
    sim.add_votable_fields("otype", "morphtype", "dimensions", "velocity", "mesDistance")
    return sim

def _make_simbad_photometry():
    sim = Simbad()
    sim.add_votable_fields("V", "B")
    return sim

SIM_BULK        = _make_simbad_bulk()
SIM_DETAIL      = _make_simbad_detail()
SIM_PHOTOMETRY  = _make_simbad_photometry()

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
# SIMBAD-STATISTIK: Einzelobjekt ausgeben (DSO-Variante)
# ============================================================ #
def _print_simbad_row(name, data):
    """Gibt die Simbad-Statistik für ein Deep-Sky-Objekt aus.
    `data` ist ein {spalte: wert}-Dict, ggf. aus mehreren Abfragen
    (Basis + Detailfelder + Photometrie) zusammengeführt."""

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

    morph = val("morph_type")
    if morph:
        lines.append(f"  Morphologie      : {morph}")

    vmag = val("V")
    bmag = val("B")
    if vmag: lines.append(f"  Ges.-Helligk.(V) : {float(vmag):.2f} mag")
    if bmag: lines.append(f"  Ges.-Helligk.(B) : {float(bmag):.2f} mag")

    maj = val("galdim_majaxis")
    min_ = val("galdim_minaxis")
    ang = val("galdim_angle")
    if maj:
        size_str = f"{float(maj):.2f}'"
        if min_ and min_ != maj:
            size_str += f" × {float(min_):.2f}'"
        if ang:
            size_str += f"  (PA {float(ang):.0f}°)"
        lines.append(f"  Winkelausdehnung : {size_str}")

    rv = val("rvz_radvel")
    if rv:
        rv_f = float(rv)
        richtung = "auf uns zu (blauverschoben)" if rv_f < 0 else "von uns weg (rotverschoben)"
        lines.append(f"  Radialgeschw.    : {rv_f:+.1f} km/s ({richtung})")

    z = val("rvz_redshift")
    if z:
        z_f = float(z)
        lines.append(f"  Rotverschiebung  : z = {z_f:.6f}")
        # Hubble-Distanz NUR bei positiver (kosmologischer) Rotverschiebung
        # sinnvoll. Negative Werte (Blauverschiebung) treten bei nahen
        # Objekten auf, deren Eigenbewegung (z.B. Lokale Gruppe, Gravitation)
        # die kosmische Expansion überwiegt (z.B. M31 nähert sich uns trotz
        # Expansion) - dort würde die Formel unsinnige negative "Distanzen"
        # liefern. Grober Näherungswert, unkorrigiert für Pekuliargeschw.
        if z_f > 0.003:
            dist_mpc = (z_f * C_KM_S) / H0_KM_S_MPC
            dist_ly = dist_mpc * 3.26156e6
            lines.append(f"  Hubble-Distanz   : ≈ {dist_mpc:,.1f} Mpc  "
                         f"(≈ {dist_ly:,.0f} Lichtjahre, aus z via H0={H0_KM_S_MPC})")

    # Literatur-Entfernung (mesDistance), sofern per query_simbad_detail
    # nachgeladen - deutlich zuverlässiger als Hubble-Distanz, v.a. für
    # nahe Objekte (Sternhaufen, galaktische Nebel, Lokale Gruppe)
    lit_dist = val("mesdistance.dist")
    lit_unit = val("mesdistance.unit")
    lit_method = val("mesdistance.method")
    if lit_dist:
        unit_str = (lit_unit or "").strip()
        method_str = f", Methode: {lit_method.strip()}" if lit_method and lit_method.strip() else ""
        lines.append(f"  Lit.-Entfernung  : {float(lit_dist):,.3f} {unit_str}{method_str}")

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
# SIMBAD: Detailfelder (1:n-Messtabelle mesDistance) nachladen
# ============================================================ #
def query_simbad_detail(name):
    """
    Holt die 1:n-Messtabelle mesDistance (publizierte Literatur-
    Entfernungen, mehrere Methoden/Einheiten möglich) für ein einzelnes
    Objekt. Da diese Tabelle pro Objekt mehrere Zeilen liefern kann, wird
    über _first_valid_per_column() für jede Spalte einzeln der erste
    gültige Wert über alle Zeilen gesucht - so geht kein vorhandener
    Messwert durch eine ungünstige Zeilen-Kombination verloren. Eine
    Verwechslung mit einem anderen Objekt ist ausgeschlossen, da hier
    immer genau ein Name abgefragt wird.
    """
    try:
        result = SIM_DETAIL.query_object(name)
    except Exception:
        return {}
    if result is None or len(result) == 0:
        return {}
    return _first_valid_per_column(result)

# ============================================================ #
# SIMBAD: Photometrie (V/B) isoliert nachladen
# ============================================================ #
def query_simbad_photometry(name):
    """
    Fragt V/B-Helligkeit separat und isoliert ab.

    HINTERGRUND: V/B lösen bei Simbad einen Inner-Join auf die
    Fluss-Tabelle aus. Objekte OHNE katalogisierte Breitbandphotometrie
    (z.B. Emissionsnebel wie M42 oder NGC 7000) liefern dann GAR KEIN
    Ergebnis mehr - auch nicht bei Einzelabfragen. Würde man V/B mit
    anderen Feldern gemeinsam abfragen, gingen dadurch alle anderen
    (durchaus vorhandenen) Basisdaten mit verloren. Durch die Isolation
    in einer eigenen Abfrage bleibt der Rest der Statistik unangetastet -
    fehlt V/B, wird einfach nur dieser Teil leer gelassen.
    """
    import warnings
    from astroquery.exceptions import NoResultsWarning
    try:
        with warnings.catch_warnings():
            # Erwarteter Fall: Objekte ohne Breitbandphotometrie (z.B. M42,
            # NGC 7000) liefern hier bewusst kein Ergebnis - keine Warnung nötig.
            warnings.simplefilter("ignore", NoResultsWarning)
            result = SIM_PHOTOMETRY.query_object(name)
    except Exception:
        return {}
    if result is None or len(result) == 0:
        return {}
    return _row_to_dict(result[0], result.colnames)

# ============================================================ #
# SIMBAD: Einzelobjekt (Basis- + Detail- + Photometriefelder kombiniert)
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
    data.update(query_simbad_photometry(name))

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
    Löst ein Deep-Sky-Objekt auf und gibt (falls vorhanden) seine Simbad-
    Statistik aus. Wird sowohl von process_constellation() als auch von
    constellation_gesamt() verwendet, damit die Detailausgabe in beiden
    Modi identisch ist.

    Basisdaten kommen aus dem Bulk-Cache (sicher zugeordnet), Detailfelder
    (Literatur-Entfernung, Photometrie) werden pro Objekt einzeln
    nachgeladen - das läuft ohnehin bereits in einem eigenen Thread und
    kann daher nicht mit einem anderen Objekt verwechselt werden.

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
    data.update(query_simbad_photometry(name))
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
    tprint(f"\n  Frage Simbad für {len(stars)} Deep-Sky-Objekte auf einmal ab...")
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
    Löst alle bekannten DSOs der Konstellation auf (inkl. Simbad-
    Detailausgabe je Objekt), berechnet den Mittelpunkt und lädt ein
    einzelnes Übersichtsbild herunter.
    """
    tprint(f"\n  Frage Simbad für {len(stars)} Deep-Sky-Objekte auf einmal ab...")
    simbad_cache = query_simbad_bulk(stars)
    tprint(f"  -> {len(simbad_cache)}/{len(stars)} Einträge erhalten.\n")

    tprint(f"  Löse Koordinaten für {len(stars)} Objekte auf...")

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
        description="Lädt Survey-Bilder für Deep-Sky-Objekte (Galaxien, Nebel, "
                    "Sternhaufen) herunter und zeigt Simbad-Statistiken."
    )
    parser.add_argument(
        "object", nargs="*",
        help="DSO-Name, z.B. M42, M31, NGC 7000, IC 434 (entfällt bei --constellation)"
    )
    parser.add_argument(
        "--fov", type=float, default=DEFAULT_FOV,
        help=f"Sichtfeld in Grad (Default: {DEFAULT_FOV}). Für kompakte Objekte "
             f"(z.B. planetarische Nebel) ggf. kleiner wählen, für ausgedehnte "
             f"Nebel/Galaxien (z.B. M31, NGC 7000) ggf. größer.",
    )
    parser.add_argument(
        "--constellation", metavar="NAME", nargs="+",
        help="Lädt alle bekannten Deep-Sky-Objekte einer Konstellation, z.B. Orion "
             "oder Ursa Major (mehrere Wörter auch ohne Anführungszeichen möglich). "
             f"Verfügbar: {', '.join(sorted(set(CONSTELLATION_DSOS)))}",
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
        stars = CONSTELLATION_DSOS.get(key)
        if stars is None:
            print(f"FEHLER: '{constellation_name}' nicht in der Liste.")
            print(f"Verfügbar: {', '.join(sorted(set(CONSTELLATION_DSOS)))}")
            sys.exit(1)
        if not stars:
            print(f"FEHLER: Für '{constellation_name}' sind aktuell keine "
                  f"bekannten Deep-Sky-Objekte hinterlegt.")
            sys.exit(1)

        formats = resolve_formats(("png",))
        print("=" * 60)
        print(f" Konstellation: {constellation_name.title()}  "
              f"({len(stars)} Deep-Sky-Objekte, Sichtfeld = {args.fov}°)")
        print("=" * 60)

        if args.gesamt:
            constellation_gesamt(constellation_name, stars, args.fov, formats, surveys, projection=args.projection)
        else:
            print(f" {MAX_WORKERS} parallele Threads")
            successes = process_constellation(stars, args.fov, formats, surveys, projection=args.projection)
            print("=" * 60)
            print(f" Fertig. {successes}/{len(stars)} Objekte erfolgreich heruntergeladen.")

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
