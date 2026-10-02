#!/usr/bin/env python3

"""
DRAUFSICHT: ERDBAHN + APOPHIS-BAHN
==================================

Blick von Norden auf die Ekliptik (heliozentrisch, J2000):
Erdbahn als Ellipse, Apophis-Bahn als echte Flugbahn
2020–2036 – mit sichtbarem Knick durch die enge
Erdnähe am 13.04.2029 (danach gilt eine neue Ellipse).
Dazu die aktuellen Positionen als Punkte. Zeitraffer
einstellbar.

Daten: Erde aus de440.bsp, Apophis aus apophis_2020_2036.bsp
(beide neben diesem Skript) – läuft offline.

Start:
    ~/astro/.venv/bin/python apophis_erde_draufsicht.py
    ~/astro/.venv/bin/python apophis_erde_draufsicht.py --speed 20
"""

import argparse
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
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
        description="Draufsicht Erdbahn + Apophis-Bahn (Ekliptik).",
    )
    p.add_argument(
        "--speed", type=float, default=5.0,
        help="Zeitraffer in Tagen pro Sekunde (Standard: 5).",
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
# Startzeit
# ============================================================

start_datetime = datetime.now(timezone.utc)

start_et = spice.str2et(
    start_datetime.strftime("%Y-%m-%d %H:%M:%S UTC")
)

start_monotonic = time.monotonic()


# ============================================================
# Bahnen abtasten (heliozentrisch, Ekliptik-Ebene)
# ============================================================

def helio_xy(target, et_liste):
    """X/Y-Positionen in AE, Sonne im Ursprung."""
    punkte = []
    for et in et_liste:
        zustand, _ = spice.spkezr(
            target, et, "ECLIPJ2000", "NONE", "SUN",
        )
        punkte.append(np.array(zustand[:3]) / AU_KM)
    return np.array(punkte)


# Erde: ein volles Jahr ab Jahresanfang (stabile Bahn)
jahr = start_datetime.year
erde_t0 = spice.str2et(f"{jahr}-01-01 00:00:00 UTC")
erde_zeiten = erde_t0 + np.linspace(0, 365.25 * TAG_S, 720)
erde_bahn = helio_xy("EARTH", erde_zeiten)

# Apophis: echte Flugbahn über die Kernel-Spanne –
# die Erdnähe 2029 ändert die Bahn, eine einzelne
# Ellipse wäre danach falsch
apo_t0 = spice.str2et("2020-01-01 00:00:00 UTC")
apo_t1 = spice.str2et("2035-12-01 00:00:00 UTC")
apophis_zeiten = np.linspace(apo_t0, apo_t1, 2922)
apophis_bahn = helio_xy(APOPHIS, apophis_zeiten)

# Enge Erdnähe als Markierung
begegnung_et = spice.str2et("2029-04-13 21:00:00 UTC")
begegnung_zustand, _ = spice.spkezr(
    APOPHIS, begegnung_et, "ECLIPJ2000", "NONE", "SUN",
)
begegnung_xy = np.array(begegnung_zustand[:3]) / AU_KM

print()
print(f"Erdbahn: {len(erde_bahn)} Punkte, "
      f"r {np.linalg.norm(erde_bahn, axis=1).min():.3f}–"
      f"{np.linalg.norm(erde_bahn, axis=1).max():.3f} AE")
print(f"Apophis-Bahn: {len(apophis_bahn)} Punkte, "
      f"r {np.linalg.norm(apophis_bahn, axis=1).min():.3f}–"
      f"{np.linalg.norm(apophis_bahn, axis=1).max():.3f} AE")


# ============================================================
# Aktuelle Positionen
# ============================================================

def aktuelle_positionen(et):
    zustand_erde, _ = spice.spkezr(
        "EARTH", et, "ECLIPJ2000", "NONE", "SUN",
    )
    zustand_apo, _ = spice.spkezr(
        APOPHIS, et, "J2000", "LT+S", "EARTH",
    )
    erde = np.array(zustand_erde[:3]) / AU_KM
    zustand_apo_helio, _ = spice.spkezr(
        APOPHIS, et, "ECLIPJ2000", "NONE", "SUN",
    )
    apo = np.array(zustand_apo_helio[:3]) / AU_KM
    distanz = float(np.linalg.norm(zustand_apo[:3])) / AU_KM
    return erde, apo, distanz


erde_now, apo_now, distanz_now = aktuelle_positionen(start_et)

print(f"Apophis–Erde jetzt: {distanz_now:.4f} AE")


# ============================================================
# Plot
# ============================================================

plt.style.use("dark_background")

fig, ax = plt.subplots(figsize=(9, 9))

ax.plot(
    erde_bahn[:, 0], erde_bahn[:, 1],
    color="tab:blue", linewidth=1.2, label="Erdbahn",
)
ax.plot(
    apophis_bahn[:, 0], apophis_bahn[:, 1],
    color="tab:orange", linewidth=1.2, label="Apophis-Bahn",
)

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

punkt_erde, = ax.plot(
    [erde_now[0]], [erde_now[1]],
    color="tab:blue", marker="o", markersize=10,
    linestyle="None", label="Erde jetzt",
)
punkt_apo, = ax.plot(
    [apo_now[0]], [apo_now[1]],
    color="tab:orange", marker="o", markersize=10,
    linestyle="None", label="Apophis jetzt",
)

ax.set_aspect("equal")
ax.set_xlabel("X [AE] (J2000-Ekliptik, heliozentrisch)")
ax.set_ylabel("Y [AE]")
ax.grid(True, alpha=0.3)
ax.legend(loc="upper right")

titel = ax.set_title(
    f"Erde + (99942) Apophis – Draufsicht – "
    f"{start_datetime.strftime('%Y-%m-%d %H:%M UTC')}  "
    f"({ARGS.speed:g} Tage/s)"
)

abstand_text = fig.text(
    0.5, 0.02,
    f"Apophis – Erde: {distanz_now * AU_KM / 1e6:8.2f} Mio km  "
    f"({distanz_now:6.4f} AE)",
    ha="center", va="bottom",
    fontsize=12, family="monospace",
)


# ============================================================
# Animation
# ============================================================

def update(frame):

    elapsed_real = time.monotonic() - start_monotonic

    current_et = start_et + elapsed_real * ARGS.speed * TAG_S

    current_datetime = start_datetime + timedelta(
        seconds=elapsed_real * ARGS.speed * TAG_S
    )

    erde, apo, distanz = aktuelle_positionen(current_et)

    punkt_erde.set_data([erde[0]], [erde[1]])
    punkt_apo.set_data([apo[0]], [apo[1]])

    titel.set_text(
        f"Erde + (99942) Apophis – Draufsicht – "
        f"{current_datetime.strftime('%Y-%m-%d %H:%M UTC')}  "
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
