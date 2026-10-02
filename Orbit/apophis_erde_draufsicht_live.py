#!/usr/bin/env python3

"""
DRAUFSICHT LIVE: ERDE + APOPHIS MIT SCHWEIF-SPUR
================================================

Wie apophis_erde_draufsicht.py, aber ohne vorgezeichnete
Bahnen: Erde und Apophis ziehen eine Schweif-Spur hinter
sich her, die nach hinten verblasst. Dadurch wird die
Bahnänderung durch die Erdnähe 2029 live sichtbar.

Daten: Erde aus de440.bsp, Apophis aus apophis_2020_2036.bsp
(beide neben diesem Skript) – läuft offline.

Am Kernel-Ende beginnt die Zeit von vorn (Spuren reset).

Start:
    ~/astro/.venv/bin/python apophis_erde_draufsicht_live.py
    ~/astro/.venv/bin/python apophis_erde_draufsicht_live.py --speed 60
"""

import argparse
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.collections import LineCollection
from datetime import datetime, timezone, timedelta
from pathlib import Path
import spiceypy as spice


# ============================================================
# Dateien
# ============================================================

BASE = Path(__file__).resolve().parent

LSK = BASE / "naif0012.tls"
SPK = BASE / "de440.bsp"
APO = BASE / "apophis_2020_2036.bsp"

# NAIF-ID (99942) Apophis (als Text, so will es spkezr)
APOPHIS = "2099942"

AU_KM = 149597870.7
TAG_S = 86400.0


# ============================================================
# Optionen
# ============================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Draufsicht live: Schweif-Spuren statt Bahnen.",
    )
    p.add_argument(
        "--speed", type=float, default=30.0,
        help="Zeitraffer in Tagen pro Sekunde (Standard: 30).",
    )
    p.add_argument(
        "--trail", type=int, default=300,
        help="Spurlänge in Punkten (Standard: 300).",
    )
    return p.parse_args()


ARGS = parse_args()


# ============================================================
# SPICE Kernel laden
# ============================================================

print("LSK:", LSK)
print("LSK vorhanden:", LSK.exists())

print("SPK:", SPK)
print("SPK vorhanden:", SPK.exists())

print("APO:", APO)
print("APO vorhanden:", APO.exists())

spice.furnsh(str(LSK))
spice.furnsh(str(SPK))
spice.furnsh(str(APO))


# ============================================================
# Kontrolle
# ============================================================

print()
print("Geladene Kernel:", spice.ktotal("ALL"))

for i in range(spice.ktotal("ALL")):
    print(spice.kdata(i, "ALL"))


# ============================================================
# Startzeit + Kernel-Ende (für den Umlauf)
# ============================================================

start_datetime = datetime.now(timezone.utc)

start_et = spice.str2et(
    start_datetime.strftime("%Y-%m-%d %H:%M:%S UTC")
)

start_monotonic = time.monotonic()

# Etwas Reserve vor Abdeckungsende (Kernel bis 12/2035)
ende_et = spice.str2et("2035-11-01 00:00:00 UTC")


# ============================================================
# Positionen
# ============================================================

def position_xy(target, et):
    zustand, _ = spice.spkezr(
        target, et, "ECLIPJ2000", "NONE", "SUN",
    )
    return np.array(zustand[:2]) / AU_KM


def distanz_au(et):
    zustand, _ = spice.spkezr(
        APOPHIS, et, "J2000", "LT+S", "EARTH",
    )
    return float(np.linalg.norm(zustand[:3])) / AU_KM


# Enge Erdnähe als feste Markierung
begegnung_et = spice.str2et("2029-04-13 21:00:00 UTC")
begegnung_xy = position_xy(APOPHIS, begegnung_et)


# ============================================================
# Plot (keine Bahnen, nur Punkte + Spuren)
# ============================================================

plt.style.use("dark_background")

fig, ax = plt.subplots(figsize=(9, 9))

ax.plot(0, 0, marker="o", color="yellow", markersize=10,
        linestyle="None", label="Sonne")

ax.plot(
    [begegnung_xy[0]], [begegnung_xy[1]],
    marker="x", color="tab:red", markersize=10,
    linestyle="None", label="Ernähe 13.04.2029",
)
ax.text(
    begegnung_xy[0], begegnung_xy[1],
    "  13.04.2029",
    color="tab:red", fontsize=9,
)

spur_erde = LineCollection([], colors="tab:blue", linewidths=1.5)
spur_apo = LineCollection([], colors="tab:orange", linewidths=1.5)
ax.add_collection(spur_erde)
ax.add_collection(spur_apo)

punkt_erde, = ax.plot(
    [], [],
    color="tab:blue", marker="o", markersize=10,
    linestyle="None", label="Erde jetzt",
)
punkt_apo, = ax.plot(
    [], [],
    color="tab:orange", marker="o", markersize=10,
    linestyle="None", label="Apophis jetzt",
)

ax.set_aspect("equal")
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
ax.set_xlabel("X [AE] (J2000-Ekliptik, heliozentrisch)")
ax.set_ylabel("Y [AE]")
ax.grid(True, alpha=0.3)
ax.legend(loc="upper right")

titel = ax.set_title("")

abstand_text = fig.text(
    0.5, 0.02,
    "",
    ha="center", va="bottom",
    fontsize=12, family="monospace",
)


# ============================================================
# Schweif-Spur mit Verlauf
# ============================================================

verlauf_erde = []
verlauf_apo = []


def spur_setzen(sammlung, verlauf, farbe):
    """Strecke aus Verlauf, nach hinten verblassend."""
    n = len(verlauf)
    if n < 2:
        sammlung.set_segments([])
        return
    punkte = np.array(verlauf[-ARGS.trail:])
    m = len(punkte)
    segmente = np.stack([punkte[:-1], punkte[1:]], axis=1)
    alpha = np.linspace(0.05, 1.0, m - 1)
    sammlung.set_segments(segmente)
    sammlung.set_colors([
        (*farbe, a) for a in alpha
    ])


FARBE_ERDE = (0.1216, 0.4667, 0.7059)   # tab:blue
FARBE_APO = (1.0, 0.4980, 0.0549)       # tab:orange


# ============================================================
# Animation
# ============================================================

sim_et = start_et
letzt_monotonic = time.monotonic()


def update(frame):

    global sim_et, letzt_monotonic

    # --------------------------------------------------------
    # Simulationszeit (mit Umlauf am Kernel-Ende)
    # --------------------------------------------------------

    jetzt = time.monotonic()
    sim_et += (jetzt - letzt_monotonic) * ARGS.speed * TAG_S
    letzt_monotonic = jetzt

    if sim_et >= ende_et:
        print("Kernel-Ende erreicht – Zeit läuft von vorn.")
        sim_et = start_et
        verlauf_erde.clear()
        verlauf_apo.clear()

    sim_datum = (
        start_datetime
        + timedelta(seconds=(sim_et - start_et))
    )

    # --------------------------------------------------------
    # Neue Positionen anhängen
    # --------------------------------------------------------

    erde_xy = position_xy("EARTH", sim_et)
    apo_xy = position_xy(APOPHIS, sim_et)

    verlauf_erde.append(erde_xy)
    verlauf_apo.append(apo_xy)

    if len(verlauf_erde) > ARGS.trail:
        del verlauf_erde[:-ARGS.trail]
    if len(verlauf_apo) > ARGS.trail:
        del verlauf_apo[:-ARGS.trail]

    distanz = distanz_au(sim_et)

    # --------------------------------------------------------
    # Zeichnen
    # --------------------------------------------------------

    spur_setzen(spur_erde, verlauf_erde, FARBE_ERDE)
    spur_setzen(spur_apo, verlauf_apo, FARBE_APO)

    punkt_erde.set_data([erde_xy[0]], [erde_xy[1]])
    punkt_apo.set_data([apo_xy[0]], [apo_xy[1]])

    titel.set_text(
        f"Erde + (99942) Apophis – live – "
        f"{sim_datum.strftime('%Y-%m-%d %H:%M UTC')}  "
        f"({ARGS.speed:g} Tage/s)"
    )

    abstand_text.set_text(
        f"Apophis – Erde: {distanz * AU_KM / 1e6:8.2f} Mio km  "
        f"({distanz:6.4f} AE)"
    )


animation = FuncAnimation(
    fig,
    update,
    interval=200,
    blit=False,
    cache_frame_data=False
)


# ============================================================
# Anzeigen
# ============================================================

plt.show()


# ============================================================
# SPICE freigeben
# ============================================================

spice.kclear()
