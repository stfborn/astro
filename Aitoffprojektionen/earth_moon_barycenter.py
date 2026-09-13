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
PCK = BASE / "pck00011.tpc"
EARTH_PCK = BASE / "earth_000101_261203_260906.bpc"


# ============================================================
# SPICE Kernel laden
# ============================================================

print("LSK:", LSK)
print("LSK vorhanden:", LSK.exists())

print("SPK:", SPK)
print("SPK vorhanden:", SPK.exists())

spice.furnsh(str(LSK))
spice.furnsh(str(SPK))
spice.furnsh(str(PCK))
spice.furnsh(str(EARTH_PCK))


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
# Erde-Mond-System
# ============================================================

EMB_DICT = {
    "Erde": "EARTH",
    "Mond": "MOON",
}


# ============================================================
# Farben
# ============================================================

BODY_COLOR_ARRAY = [
    "tab:blue",
    "lightgray",
]


# ============================================================
# Positionen bestimmen
# ============================================================

def calculate_positions(et):

    positions = {}

    # --------------------------------------------------------
    # Erde und Mond relativ zum Erde-Mond-Baryzentrum
    # --------------------------------------------------------

    for body_name, body in EMB_DICT.items():

        state, _ = spice.spkezr(
            body,
            et,
            "J2000",
            "NONE",
            "EARTH BARYCENTER"
        )

        position = state[:3]
        velocity = state[3:6]

        distance = np.linalg.norm(position)
        velocity_magnitude = np.linalg.norm(velocity)

        # ----------------------------------------------------
        # Kartesische Koordinaten
        # ----------------------------------------------------

        x, y, z = position

        # ----------------------------------------------------
        # Rektaszension
        # ----------------------------------------------------

        ra = np.arctan2(y, x)

        ra = (
            (ra + np.pi) % (2 * np.pi)
            - np.pi
        )

        # ----------------------------------------------------
        # Deklination
        # ----------------------------------------------------

        dec = np.arcsin(
            z / distance
        )

        # ----------------------------------------------------
        # Werte speichern
        # ----------------------------------------------------

        positions[body_name] = {
            "ra": ra,
            "dec": dec,
            "distance": distance,
            "velocity": velocity_magnitude,
        }

    # --------------------------------------------------------
    # Erde-Mond-Distanz
    # --------------------------------------------------------

    earth_state, _ = spice.spkezr(
        "EARTH",
        et,
        "J2000",
        "NONE",
        "EARTH BARYCENTER"
    )

    moon_state, _ = spice.spkezr(
        "MOON",
        et,
        "J2000",
        "NONE",
        "EARTH BARYCENTER"
    )

    earth_moon_vector = (
        moon_state[:3]
        - earth_state[:3]
    )

    earth_moon_distance = np.linalg.norm(
        earth_moon_vector
    )

    positions["Erde-Mond-Distanz"] = (
        earth_moon_distance
    )

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
    f"Erde-Mond-System – Bezugspunkt "
    f"Erde-Mond-Baryzentrum – "
    f"{start_datetime.strftime('%Y-%m-%d %H:%M:%S UTC')}",
    fontsize=12,
    fontweight="bold",
    x=0.50,
    y=1.06
)


# ============================================================
# Erde und Mond zeichnen
# ============================================================

body_points = {}

for body_name, body_color in zip(
    EMB_DICT,
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
        label=body_name
    )

    body_points[body_name] = point


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
# Infobox
# ============================================================

info_lines = []

for body_name in EMB_DICT:

    data = positions[body_name]

    ra_deg = np.degrees(
        data["ra"]
    )

    dec_deg = np.degrees(
        data["dec"]
    )

    info_lines.append(
        f"{body_name:5s} "
        f"Baryzentrum {data['distance']:9.2f} km   "
        f"v {data['velocity']:6.3f} km/s   "
        f"RA {ra_deg:7.2f}°   "
        f"Dec {dec_deg:6.2f}°"
    )


info_lines.append(
    f"Erde-Mond-Distanz "
    f"{positions['Erde-Mond-Distanz']:9.2f} km"
)


info_text = fig.text(
    0.5,
    0.055,
    "\n".join(info_lines),
    transform=fig.transFigure,
    ha="center",
    va="bottom",
    fontsize=11,
    family="monospace"
)


# ============================================================
# Gitter
# ============================================================

plt.grid(True)


# ============================================================
# Ekliptik
# ============================================================

ra_ecliptic = np.linspace(
    -np.pi,
    np.pi,
    1000
)

# Schiefe der Ekliptik im J2000-System

epsilon = np.radians(
    23.4392911
)

dec_ecliptic = np.arcsin(
    np.sin(epsilon)
    * np.sin(ra_ecliptic)
)

ax.plot(
    ra_ecliptic,
    dec_ecliptic,
    linestyle="--",
    linewidth=1.2,
    color="white",
    alpha=0.7,
    label="Ekliptik"
)


# ============================================================
# Legende
# ============================================================

ax.legend(
    loc="upper right",
    bbox_to_anchor=(1.15, 1.02)
)


# ============================================================
# Scaling / Layout
# ============================================================

plt.subplots_adjust(
    bottom=0.25,
    top=0.92
)


# ============================================================
# Animation
# ============================================================

def update(frame):

    # --------------------------------------------------------
    # Tatsächlich vergangene Zeit
    # --------------------------------------------------------

    elapsed_seconds = (
        time.monotonic()
        - start_monotonic
    )

    # --------------------------------------------------------
    # 1 reale Sekunde = 1 Sekunde Simulation
    # --------------------------------------------------------

    current_et = (
        start_et
        + elapsed_seconds
    )

    current_datetime = (
        start_datetime
        + timedelta(seconds=elapsed_seconds)
    )


    # --------------------------------------------------------
    # Neue Positionen
    # --------------------------------------------------------

    positions = calculate_positions(
        current_et
    )


    # --------------------------------------------------------
    # Erde und Mond aktualisieren
    # --------------------------------------------------------

    for body_name in EMB_DICT:

        ra = positions[body_name]["ra"]
        dec = positions[body_name]["dec"]

        body_points[body_name].set_data(
            [ra],
            [dec]
        )


    # --------------------------------------------------------
    # Titel aktualisieren
    # --------------------------------------------------------

    title.set_text(
        f"Erde-Mond-System – Bezugspunkt "
        f"Erde-Mond-Baryzentrum – "
        f"{current_datetime.strftime('%Y-%m-%d %H:%M:%S UTC')}"
    )


    # --------------------------------------------------------
    # Infobox aktualisieren
    # --------------------------------------------------------

    info_lines = []

    for body_name in EMB_DICT:

        data = positions[body_name]

        ra_deg = np.degrees(
            data["ra"]
        )

        dec_deg = np.degrees(
            data["dec"]
        )

        info_lines.append(
            f"{body_name:5s} "
            f"Baryzentrum {data['distance']:9.2f} km   "
            f"v {data['velocity']:6.3f} km/s   "
            f"RA {ra_deg:7.2f}°   "
            f"Dec {dec_deg:6.2f}°"
        )


    info_lines.append(
        f"Erde-Mond-Distanz "
        f"{positions['Erde-Mond-Distanz']:9.2f} km"
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
