#!/usr/bin/env python3

"""
JUPITER ÜBER MAINZ
===============

Skyfield + JPL DE421

Berechnet für Mainz (Standard: 49.9858° N, 8.2791° E),
wann Jupiter wo steht:

- Aufgang / Kulmination (höchster Stand) / Untergang
- Höhe / Azimut / Himmelsrichtung im Tagesverlauf
- Rektaszension / Deklination / Entfernung
- aktuelle Position ("wo steht Jupiter jetzt?")

Start:
    python3 jupiter_mainz.py --help
    python3 jupiter_mainz.py
    python3 jupiter_mainz.py --datum 2026-09-28 --tage 7
    python3 jupiter_mainz.py --verlauf 2026-09-28 --schritt 60
"""

import argparse
from datetime import datetime, timedelta
from pathlib import Path
from zoneinfo import ZoneInfo

from skyfield.api import load, wgs84


# ============================================================
# STANDORT MAINZ
# ============================================================

MAINZ_LAT = 49.9858
MAINZ_LON = 8.2791
MAINZ_ELEVATION = 90

# Schwelle für enge Mond–Jupiter-Konjunktion (Grad)
MOND_SCHWELLE = 5.0

BERLIN = ZoneInfo("Europe/Berlin")

BASE = Path(__file__).resolve().parent
STANDARD_EPHEMERIDE = BASE / "de421.bsp"


# ============================================================
# HILFSFUNKTIONEN
# ============================================================

def himmelsrichtung(azimut_grad):
    """Azimut 0°=N ... 360°=N -> 16-teilige Kompassrichtung."""
    namen = [
        "N", "NNO", "NO", "ONO",
        "O", "OSO", "SO", "SSO",
        "S", "SSW", "SW", "WSW",
        "W", "WNW", "NW", "NNW",
    ]
    index = int((azimut_grad + 11.25) // 22.5) % 16
    return namen[index]


def fmt_zeit(dt_lokal):
    return dt_lokal.strftime("%H:%M")


def fmt_datum(dt_lokal):
    return dt_lokal.strftime("%d.%m.%Y")


# ============================================================
# POSITION
# ============================================================

def jupiter_position(mainz, jupiter, t):
    """Höhe/Azimut/Distanz von Mainz aus, dazu RA/Dec."""
    scheinbar = mainz.at(t).observe(jupiter).apparent()
    hoehe, azimut, distanz = scheinbar.altaz()
    ra, dec, _ = scheinbar.radec()
    return {
        "hoehe": hoehe.degrees,
        "azimut": azimut.degrees % 360.0,
        "distanz_au": distanz.au,
        "ra_h": ra.hours,
        "dec_deg": dec.degrees,
    }


# ============================================================
# EIN TAG: AUFGANG / KULMINATION / UNTERGANG
# ============================================================

def tagesereignisse(ts, mainz, jupiter, lokales_datum):
    """
    Sucht im lokalen Kalendertag (00:00-24:00 Europe/Berlin)
    nach Aufgang, Kulmination und Untergang.

    Abtastung: 1 Minute. Auf-/Untergang = Nulldurchgang
    der Höhe, Kulmination = Maximum der Höhe.
    """
    start_lokal = datetime(
        lokales_datum.year, lokales_datum.month, lokales_datum.day,
        0, 0, tzinfo=BERLIN,
    )
    minuten = [start_lokal + timedelta(minutes=m) for m in range(24 * 60 + 1)]
    zeiten = ts.from_datetimes(minuten)

    hoehen = []
    azimute = []
    for t in zeiten:
        pos = jupiter_position(mainz, jupiter, t)
        hoehen.append(pos["hoehe"])
        azimute.append(pos["azimut"])

    aufgang = None
    untergang = None

    for i in range(1, len(hoehen)):
        vorher, nachher = hoehen[i - 1], hoehen[i]
        # Aufgang: von unter auf über Horizont
        if vorher <= 0 < nachher:
            # linear interpolieren für genauere Zeit
            bruch = (0 - vorher) / (nachher - vorher)
            aufgang = minuten[i - 1] + timedelta(minutes=bruch)
        # Untergang: von über auf unter Horizont
        if vorher > 0 >= nachher:
            bruch = (vorher - 0) / (vorher - nachher)
            untergang = minuten[i - 1] + timedelta(minutes=bruch)

    # Kulmination = höchste Höhe des Tages
    k_index = max(range(len(hoehen)), key=lambda i: hoehen[i])
    kulmination = minuten[k_index]
    k_pos = jupiter_position(mainz, jupiter, zeiten[k_index])

    return {
        "datum": lokales_datum,
        "aufgang": aufgang,
        "kulmination": kulmination,
        "kulmi_hoehe": k_pos["hoehe"],
        "kulmi_azimut": k_pos["azimut"],
        "kulmi_ra": k_pos["ra_h"],
        "kulmi_dec": k_pos["dec_deg"],
        "kulmi_dist": k_pos["distanz_au"],
        "untergang": untergang,
    }


# ============================================================
# TAGESVERLAUF
# ============================================================

def tagesverlauf(ts, mainz, jupiter, lokales_datum, schritt_minuten):
    """Liste von Positionen über den Tag im Schritt-Raster."""
    start_lokal = datetime(
        lokales_datum.year, lokales_datum.month, lokales_datum.day,
        0, 0, tzinfo=BERLIN,
    )
    ergebnis = []
    m = 0
    while m <= 24 * 60:
        dt_lokal = start_lokal + timedelta(minutes=m)
        t = ts.from_datetime(dt_lokal)
        pos = jupiter_position(mainz, jupiter, t)
        ergebnis.append((dt_lokal, pos))
        m += schritt_minuten
    return ergebnis


# ============================================================
# ENGE MOND-KONJUNKTION
# ============================================================

def mond_konjunktionen(ts, erde, mond, jupiter, start_datum, tage, schwelle):
    """
    Sucht lokale Minima des Winkelabstands Mond–Jupiter
    unterhalb der Schwelle (Abtastung: 1 Stunde).
    """
    zeiten_lokal = []
    t = datetime(
        start_datum.year, start_datum.month, start_datum.day,
        0, 0, tzinfo=BERLIN,
    )
    ende = t + timedelta(days=tage)
    while t <= ende:
        zeiten_lokal.append(t)
        t += timedelta(hours=1)

    sky_zeiten = ts.from_datetimes(zeiten_lokal)
    abstaende = []
    for st in sky_zeiten:
        mond_pos = erde.at(st).observe(mond).apparent()
        jupiter_pos = erde.at(st).observe(jupiter).apparent()
        abstaende.append(mond_pos.separation_from(jupiter_pos).degrees)

    ereignisse = []
    for i in range(1, len(abstaende) - 1):
        if abstaende[i] <= abstaende[i - 1] and abstaende[i] <= abstaende[i + 1]:
            if abstaende[i] <= schwelle:
                ereignisse.append((zeiten_lokal[i], abstaende[i]))
    return ereignisse


# ============================================================
# AUSGABE
# ============================================================

def print_uebersicht(ereignisse):
    print()
    print("=" * 78)
    print("JUPITER ÜBER MAINZ – AUFGANG / KULMINATION / UNTERGANG")
    print("=" * 78)
    print()
    print(
        f"{'Datum':<12} {'Aufgang':<8} {'Kulmination':<22} {'Untergang':<8}"
    )
    print("-" * 78)
    for e in ereignisse:
        auf = fmt_zeit(e["aufgang"]) if e["aufgang"] else "--:--"
        unt = fmt_zeit(e["untergang"]) if e["untergang"] else "--:--"
        kul = (
            f"{fmt_zeit(e['kulmination'])} "
            f"({e['kulmi_hoehe']:4.1f}°, "
            f"{e['kulmi_azimut']:5.1f}°/{himmelsrichtung(e['kulmi_azimut'])})"
        )
        print(f"{fmt_datum(e['datum']):<12} {auf:<8} {kul:<22} {unt:<8}")
    print()
    print("Höhe = über Horizont bei Kulmination, Azimut 0°=N 90°=O 180°=S 270°=W.")
    print("Zeiten in Europe/Berlin (MESZ/MEZ).")


def print_verlauf(verlauf):
    print()
    print("=" * 78)
    print(
        f"JUPITER-VERLAUF – {fmt_datum(verlauf[0][0])} "
        f"(Europe/Berlin, alle Zeiten lokal)"
    )
    print("=" * 78)
    print()
    print(
        f"{'Uhr':<6} {'Höhe':<8} {'Azimut':<14} "
        f"{'RA':<9} {'Dec':<9} {'Entf.':<9} Status"
    )
    print("-" * 78)
    for dt_lokal, pos in verlauf:
        status = "sichtbar" if pos["hoehe"] > 0 else "unt. Horizont"
        print(
            f"{fmt_zeit(dt_lokal):<6} "
            f"{pos['hoehe']:6.1f}° "
            f"{pos['azimut']:6.1f}°/{himmelsrichtung(pos['azimut']):<3} "
            f"{pos['ra_h']:5.2f}h "
            f"{pos['dec_deg']:6.1f}° "
            f"{pos['distanz_au']:5.3f}AE "
            f"{status}"
        )
    print()


def print_konjunktionen(ereignisse):
    print()
    print("=" * 78)
    print("ENGE MOND–JUPITER-KONJUNKTIONEN")
    print("=" * 78)
    print()
    if not ereignisse:
        print("Keine enge Konjunktion im Zeitraum.")
    for zeit, abstand in ereignisse:
        print(
            f"\033[31m  {zeit.strftime('%d.%m.%Y %H:%M')}  "
            f"Mond–Jupiter {abstand:5.2f}°\033[0m"
        )
    print()


def print_jetzt(ts, mainz, jupiter):
    jetzt_lokal = datetime.now(BERLIN)
    t = ts.from_datetime(jetzt_lokal)
    pos = jupiter_position(mainz, jupiter, t)
    print()
    print("=" * 78)
    print(
        f"JUPITER JETZT – {jetzt_lokal.strftime('%d.%m.%Y %H:%M')} "
        f"Europe/Berlin (Mainz)"
    )
    print("=" * 78)
    print()
    print(f"  Höhe:               {pos['hoehe']:6.1f}°")
    print(
        f"  Azimut:             {pos['azimut']:6.1f}° "
        f"({himmelsrichtung(pos['azimut'])})"
    )
    print(f"  Rektaszension:      {pos['ra_h']:6.2f} h")
    print(f"  Deklination:        {pos['dec_deg']:6.1f}°")
    print(f"  Entfernung:         {pos['distanz_au']:6.3f} AE "
          f"({pos['distanz_au'] * 149.5978707:6.1f} Mio km)")
    if pos["hoehe"] > 0:
        print("  Status:             über Horizont, sichtbar")
    else:
        print("  Status:             unter Horizont")
    print()


# ============================================================
# HAUPTPROGRAMM
# ============================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Jupiter-Stand für Mainz: Aufgang/Kulmination/Untergang "
                    "und Tagesverlauf (Skyfield + DE421).",
    )
    p.add_argument(
        "--datum", default=None,
        help="Startdatum JJJJ-MM-TT (Standard: heute, Europe/Berlin).",
    )
    p.add_argument(
        "--tage", type=int, default=7,
        help="Anzahl Tage ab Startdatum für die Übersicht (Standard: 7).",
    )
    p.add_argument(
        "--schritt", type=int, default=60,
        help="Schrittweite in Minuten für --verlauf (Standard: 60).",
    )
    p.add_argument(
        "--verlauf", default=None, metavar="JJJJ-MM-TT",
        help="Ausführlicher Tagesverlauf für dieses Datum "
             "(Standard: Startdatum). Mit --kein-verlauf abschalten.",
    )
    p.add_argument(
        "--kein-verlauf", action="store_true",
        help="Nur Übersicht + Jetzt, kein Stundentableau.",
    )
    p.add_argument("--lat", type=float, default=MAINZ_LAT)
    p.add_argument("--lon", type=float, default=MAINZ_LON)
    p.add_argument("--elev", type=float, default=MAINZ_ELEVATION)
    p.add_argument(
        "--ephemeride", default=str(STANDARD_EPHEMERIDE),
        help="Pfad zur DE421-BSP-Datei.",
    )
    return p.parse_args()


def main():
    args = parse_args()

    if args.datum:
        start_datum = datetime.strptime(args.datum, "%Y-%m-%d").date()
    else:
        start_datum = datetime.now(BERLIN).date()

    eph_path = Path(args.ephemeride)
    if not eph_path.exists():
        raise SystemExit(f"Ephemeride nicht gefunden: {eph_path}")

    print(f"Lade Ephemeride {eph_path} ...")
    ts = load.timescale()
    eph = load(str(eph_path))

    erde = eph["earth"]
    jupiter = eph["jupiter barycenter"]
    mond = eph["moon"]
    mainz = erde + wgs84.latlon(args.lat, args.lon, elevation_m=args.elev)

    print(f"Standort: {args.lat:.4f}° N, {args.lon:.4f}° E")

    # ---- Übersicht über N Tage ----
    ereignisse = []
    for d in range(args.tage):
        tag = start_datum + timedelta(days=d)
        ereignisse.append(tagesereignisse(ts, mainz, jupiter, tag))
    print_uebersicht(ereignisse)

    # ---- Enge Mond-Konjunktion (rot) ----
    konjunktionen = mond_konjunktionen(
        ts, erde, mond, jupiter, start_datum, args.tage, MOND_SCHWELLE,
    )
    print_konjunktionen(konjunktionen)

    # ---- Tagesverlauf ----
    if not args.kein_verlauf:
        if args.verlauf:
            verlauf_datum = datetime.strptime(args.verlauf, "%Y-%m-%d").date()
        else:
            verlauf_datum = start_datum
        verlauf = tagesverlauf(ts, mainz, jupiter, verlauf_datum, args.schritt)
        print_verlauf(verlauf)

    # ---- Jetzt ----
    print_jetzt(ts, mainz, jupiter)


if __name__ == "__main__":
    main()
