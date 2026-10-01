#!/usr/bin/env python3

"""
MOND–PLANET-LINIEN ÜBER MAINZ
=============================

Skyfield + JPL DE421

Zeichnet über einen Zeitraum den Winkelabstand
Mond–Mars, Mond–Venus, Mond–Jupiter als Linien.

Enge Begegnungen (Konjunktionen) erkennt man dort,
wo eine Linie ein Minimum hat bzw. die Schwellenlinie
(z.B. 5°) unterschreitet.

Hinweis: Echte "Oppositionen" (Sonne–Erde–Planet = 180°)
gibt es nur bei äußeren Planeten (Mars, Jupiter) – nicht
bei Venus und nicht mit dem Mond. Hier geht es um enge
Mond–Planet-Begegnungen am Himmel.

Start:
    python3 mond_planet_linien.py --help
    python3 mond_planet_linien.py --datum 2026-09-28 --tage 30
    python3 mond_planet_linien.py --tage 30 --schwelle 5 --output begegnungen.png
"""

import argparse
from datetime import datetime, timedelta
from pathlib import Path
from zoneinfo import ZoneInfo

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from skyfield.api import load


# ============================================================
# EINSTELLUNGEN
# ============================================================

BERLIN = ZoneInfo("Europe/Berlin")

BASE = Path(__file__).resolve().parent
STANDARD_EPHEMERIDE = BASE / "de421.bsp"

PLANETEN_KEYS = {
    "Mars": "mars",
    "Venus": "venus",
    "Jupiter": "jupiter barycenter",
}

FARBEN = {
    "Mars": "tab:red",
    "Venus": "tab:orange",
    "Jupiter": "tab:blue",
}


# ============================================================
# BERECHNUNG
# ============================================================

def lade_ephemeride(pfad):
    eph_path = Path(pfad)
    if not eph_path.exists():
        raise SystemExit(f"Ephemeride nicht gefunden: {eph_path}")
    print(f"Lade Ephemeride {eph_path} ...")
    ts = load.timescale()
    eph = load(str(eph_path))
    return ts, eph


def abstaende(ts, eph, start_datum, tage, schritt_stunden):
    """Winkelabstand Mond–Planet (geozentrisch) im Zeitraster."""
    erde = eph["earth"]
    mond = eph["moon"]
    koerper = {name: eph[key] for name, key in PLANETEN_KEYS.items()}

    zeiten_lokal = []
    t = datetime(
        start_datum.year, start_datum.month, start_datum.day,
        0, 0, tzinfo=BERLIN,
    )
    ende = t + timedelta(days=tage)
    while t <= ende:
        zeiten_lokal.append(t)
        t += timedelta(hours=schritt_stunden)

    sky_zeiten = ts.from_datetimes(zeiten_lokal)
    reihen = {name: [] for name in koerper}

    for st in sky_zeiten:
        mond_pos = erde.at(st).observe(mond).apparent()
        for name, body in koerper.items():
            planet_pos = erde.at(st).observe(body).apparent()
            reihen[name].append(mond_pos.separation_from(planet_pos).degrees)

    return zeiten_lokal, reihen


def begegnungen(zeiten_lokal, reihen, schwelle_grad):
    """
    Lokale Minima unterhalb der Schwelle suchen.
    Gibt Liste (planet, zeitpunkt, abstand) zurück.
    """
    ereignisse = []
    for name, werte in reihen.items():
        for i in range(1, len(werte) - 1):
            if werte[i] <= werte[i - 1] and werte[i] <= werte[i + 1]:
                if werte[i] <= schwelle_grad:
                    ereignisse.append((name, zeiten_lokal[i], werte[i]))
    return sorted(ereignisse, key=lambda e: e[1])


# ============================================================
# DIAGRAMM
# ============================================================

def zeichne(zeiten_lokal, reihen, ereignisse, schwelle_grad,
            start_datum, tage, output):
    fig, ax = plt.subplots(figsize=(12, 6))

    for name, werte in reihen.items():
        ax.plot(zeiten_lokal, werte, label=f"Mond–{name}",
                color=FARBEN.get(name))

    ax.axhline(schwelle_grad, color="gray", linestyle="--",
               linewidth=1, label=f"Schwelle {schwelle_grad}°")

    for k, (name, zeit, abstand) in enumerate(ereignisse):
        ax.plot(zeit, abstand, marker="o",
                color=FARBEN.get(name, "black"))
        # gestaffelt, damit nahe Termine sich nicht überdecken
        hoehe = 12 if k % 2 == 0 else 26
        ax.annotate(
            f"{name} {abstand:.1f}°\n{zeit.strftime('%d.%m.')}",
            xy=(zeit, abstand),
            xytext=(0, hoehe),
            textcoords="offset points",
            ha="center",
            fontsize=8,
        )

    ax.set_ylim(0, max(60.0, schwelle_grad + 10))
    ax.set_ylabel("Winkelabstand Mond–Planet [Grad]")
    ax.set_xlabel("Datum (Europe/Berlin)")
    ax.set_title(
        f"Mond–Planet-Begegnungen ab {start_datum.strftime('%d.%m.%Y')} "
        f"({tage} Tage) – Minima = enge Begegnungen"
    )
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d.%m."))
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(output, dpi=150)
    print(f"Diagramm gespeichert: {output}")


# ============================================================
# HAUPTPROGRAMM
# ============================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Linien-Diagramm Mond–Mars/Venus/Jupiter: "
                    "Minima zeigen enge Begegnungen.",
    )
    p.add_argument("--datum", default=None,
                   help="Startdatum JJJJ-MM-TT (Standard: heute).")
    p.add_argument("--tage", type=int, default=30,
                   help="Zeitraum in Tagen (Standard: 30).")
    p.add_argument("--schritt", type=float, default=3.0,
                   help="Raster in Stunden (Standard: 3).")
    p.add_argument("--schwelle", type=float, default=5.0,
                   help="Schwelle in Grad für Markierung (Standard: 5).")
    p.add_argument("--output", default=None,
                   help="PNG-Datei (Standard: mond_planet_linien_DATUM.png).")
    p.add_argument("--ephemeride", default=str(STANDARD_EPHEMERIDE))
    return p.parse_args()


def main():
    args = parse_args()

    if args.datum:
        start_datum = datetime.strptime(args.datum, "%Y-%m-%d").date()
    else:
        start_datum = datetime.now(BERLIN).date()

    output = args.output or (
        f"mond_planet_linien_{start_datum.strftime('%Y-%m-%d')}.png"
    )

    ts, eph = lade_ephemeride(args.ephemeride)
    zeiten_lokal, reihen = abstaende(
        ts, eph, start_datum, args.tage, args.schritt,
    )
    ereignisse = begegnungen(zeiten_lokal, reihen, args.schwelle)

    print()
    print("ENGE MOND–PLANET-BEGEGNUNGEN "
          f"(Schwelle {args.schwelle}°):")
    print("-" * 60)
    if not ereignisse:
        print("Keine Begegnung unter der Schwelle im Zeitraum.")
    for name, zeit, abstand in ereignisse:
        print(f"  {zeit.strftime('%d.%m.%Y %H:%M')}  "
              f"Mond–{name:<7} {abstand:5.2f}°")
    print()

    zeichne(zeiten_lokal, reihen, ereignisse, args.schwelle,
            start_datum, args.tage, output)


if __name__ == "__main__":
    main()
