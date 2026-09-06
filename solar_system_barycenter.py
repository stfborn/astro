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


# ============================================================
# SPICE Kernel laden
# ============================================================

print("LSK:", LSK)
print("LSK vorhanden:", LSK.exists())

print("SPK:", SPK)
print("SPK vorhanden:", SPK.exists())

spice.furnsh(str(LSK))
spice.furnsh(str(SPK))


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
# Planeten
# ============================================================

SOLSYS_DICT = {
    "Merkur": "MERCURY",
    "Venus": "VENUS",
    "Erde": "EARTH",
    "Mars": "MARS BARYCENTER",
    "Jupiter": "JUPITER BARYCENTER",
    "Saturn": "SATURN BARYCENTER",
    "Uranus": "URANUS BARYCENTER",
    "Neptun": "NEPTUNE BARYCENTER",
}


# ============================================================
# Farben
# ============================================================

BODY_COLOR_ARRAY = [
    "tab:gray",
    "tab:orange",
    "tab:blue",
    "tab:red",
    "tab:brown",
    "tab:pink",
    "tab:cyan",
    "tab:green",
]


# ============================================================
# Positionen bestimmen
# ============================================================

def calculate_positions(et):

    positions = {}

    for body_name, body in SOLSYS_DICT.items():

        # ----------------------------------------------------
        # Position relativ zum SSB
        # ----------------------------------------------------

        state, light_time = spice.spkezr(
            body,
            et,
            "J2000",
            "LT+S",
            "SSB"
        )

        x, y, z = state[:3]

        distance_ssb = np.linalg.norm(state[:3])

        # ----------------------------------------------------
        # Position relativ zur Sonne
        # ----------------------------------------------------

        state_sun, _ = spice.spkezr(
            body,
            et,
            "J2000",
            "LT+S",
            "SUN"
        )

        distance_sun = np.linalg.norm(state_sun[:3])

        velocity_sun = np.linalg.norm(state_sun[3:6])

        # ----------------------------------------------------
        # Position relativ zur Erde
        # ----------------------------------------------------

        state_earth, _ = spice.spkezr(
            body,
            et,
            "J2000",
            "LT+S",
            "EARTH"
        )

        distance_earth = np.linalg.norm(state_earth[:3])

        # ----------------------------------------------------
        # Rektaszension
        # ----------------------------------------------------

        ra = np.arctan2(y, x)

        ra = (ra + np.pi) % (2 * np.pi) - np.pi

        # ----------------------------------------------------
        # Deklination
        # ----------------------------------------------------

        dec = np.arcsin(z / distance_ssb)

        # ----------------------------------------------------
        # Werte speichern
        # ----------------------------------------------------

        positions[body_name] = {
            "ra": ra,
            "dec": dec,
            "distance": distance_ssb,
            "distance_sun": distance_sun,
            "distance_earth": distance_earth,
            "velocity_sun": velocity_sun,
        }

    return positions


# ============================================================
# Erste Position
# ============================================================

positions = calculate_positions(start_et)


# ============================================================
# Aitoff-Plot
# ============================================================

plt.style.use("dark_background")

fig = plt.figure(figsize=(12, 8))

ax = plt.subplot(
    projection="aitoff"
)


# ============================================================
# Titel
# ============================================================

title = plt.title(
    f"Schweremittelpunkt Sonnensystem – "
    f"{start_datetime.strftime('%Y-%m-%d %H:%M:%S UTC')}",
    fontsize=12,
    fontweight="bold",
    x=0.50,
    y=1.06
)


# ============================================================
# Planeten zeichnen
# ============================================================

planet_points = {}

for body_name, body_color in zip(
    SOLSYS_DICT,
    BODY_COLOR_ARRAY
):

    ra = positions[body_name]["ra"]
    dec = positions[body_name]["dec"]

    point, = ax.plot(
        [ra],
        [dec],
        color=body_color,
        marker="o",
        linestyle="None",
        markersize=10,
        label=body_name.capitalize()
    )

    planet_points[body_name] = point


# ============================================================
# Achsen
# ============================================================

plt.xticks(
    ticks=np.radians(
        [
            -150, -120, -90, -60, -30,
            0,
            30, 60, 90, 120, 150
        ]
    ),
    labels=[
        "10 h",
        "8 h",
        "6 h",
        "4 h",
        "2 h",
        "0 h",
        "22 h",
        "20 h",
        "18 h",
        "16 h",
        "14 h"
    ]
)

plt.xlabel("Rektaszension")
plt.ylabel("Deklination")


# ============================================================
# Legende
# ============================================================

plt.legend(
    loc="upper right",
    bbox_to_anchor=(1.15, 1.02)
)


# ============================================================
# Infobox
# ============================================================

info_lines = []

for body_name in SOLSYS_DICT:

    data = positions[body_name]

    ra_deg = np.degrees(data["ra"])
    dec_deg = np.degrees(data["dec"])

    info_lines.append(
        f"{body_name.capitalize():8s} "
        f"SSB {data['distance']/1e6:7.2f} Mio km   "
        f"Sonne {data['distance_sun']/1e6:7.2f} Mio km   "
        f"Erde {data['distance_earth']/1e6:7.2f} Mio km   "
        f"v {data['velocity_sun']:5.2f} km/s   "
        f"RA {ra_deg:7.2f}°   "
        f"Dec {dec_deg:6.2f}°"
    )


info_text = ax.text(
    0.5,
    -0.12,
    "\n".join(info_lines),
    transform=ax.transAxes,
    ha="center",
    va="top",
    fontsize=11,
    family="monospace"
)


# ============================================================
# Gitter
# ============================================================

plt.grid(True)


# ============================================================
# Scaling / Layout unverändert
# ============================================================

plt.subplots_adjust(
    bottom=0.32,
    top=0.92
)


# ============================================================
# Animation
# ============================================================

def update(frame):

    # --------------------------------------------------------
    # Tatsächlich vergangene Zeit
    # --------------------------------------------------------

    elapsed_seconds = time.monotonic() - start_monotonic

    # --------------------------------------------------------
    # 1 reale Sekunde = 1 Sekunde Simulation
    # --------------------------------------------------------

    current_et = start_et + elapsed_seconds

    current_datetime = (
        start_datetime
        + timedelta(seconds=elapsed_seconds)
    )

    # --------------------------------------------------------
    # Neue Positionen
    # --------------------------------------------------------

    positions = calculate_positions(current_et)

    # --------------------------------------------------------
    # Planeten aktualisieren
    # --------------------------------------------------------

    for body_name in SOLSYS_DICT:

        ra = positions[body_name]["ra"]
        dec = positions[body_name]["dec"]

        planet_points[body_name].set_data(
            [ra],
            [dec]
        )

    # --------------------------------------------------------
    # Titel aktualisieren
    # --------------------------------------------------------

    title.set_text(
        f"Schweremittelpunkt Sonnensystem – "
        f"{current_datetime.strftime('%Y-%m-%d %H:%M:%S UTC')}"
    )

    # --------------------------------------------------------
    # Infobox aktualisieren
    # --------------------------------------------------------

    info_lines = []

    for body_name in SOLSYS_DICT:

        data = positions[body_name]

        ra_deg = np.degrees(data["ra"])
        dec_deg = np.degrees(data["dec"])

        info_lines.append(
            f"{body_name.capitalize():8s} "
            f"SSB {data['distance']/1e6:7.2f} Mio km   "
            f"Sonne {data['distance_sun']/1e6:7.2f} Mio km   "
            f"Erde {data['distance_earth']/1e6:7.2f} Mio km   "
            f"v {data['velocity_sun']:5.2f} km/s   "
            f"RA {ra_deg:7.2f}°   "
            f"Dec {dec_deg:6.2f}°"
        )

    info_text.set_text(
        "\n".join(info_lines)
    )


# ============================================================
# Animation starten
# ============================================================

animation = FuncAnimation(
    fig,
    update,
    interval=2000,
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
