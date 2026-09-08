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

    # --------------------------------------------------------
    # Erdmittelpunkt relativ zum SSB
    # --------------------------------------------------------

    earth_state, _ = spice.spkezr(
        "EARTH",
        et,
        "J2000",
        "NONE",
        "SSB"
    )

    earth_ssb = earth_state[:3]

    # --------------------------------------------------------
    # Planeten
    # --------------------------------------------------------

    for body_name, body in SOLSYS_DICT.items():

        # ----------------------------------------------------
        # Position des Körpers relativ zum SSB
        # ----------------------------------------------------

        state, _ = spice.spkezr(
            body,
            et,
            "J2000",
            "NONE",
            "SSB"
        )

        target_ssb = state[:3]

        # ----------------------------------------------------
        # Position des Körpers relativ zum Erdmittelpunkt
        # ----------------------------------------------------

        earth_to_body = (
            target_ssb - earth_ssb
        )

        distance_earth = np.linalg.norm(
            earth_to_body
        )

        # ----------------------------------------------------
        # Erde selbst
        # Bezugspunkt = Erdmittelpunkt
        # ----------------------------------------------------

        if distance_earth == 0.0:

            state_sun, _ = spice.spkezr(
                "EARTH",
                et,
                "J2000",
                "LT+S",
                "SUN"
            )

            distance_sun = np.linalg.norm(
                state_sun[:3]
            )

            velocity_sun = np.linalg.norm(
                state_sun[3:6]
            )

            positions[body_name] = {
                "ra": np.nan,
                "dec": np.nan,
                "distance": 0.0,
                "distance_sun": distance_sun,
                "distance_earth": 0.0,
                "velocity_sun": velocity_sun,
            }

            continue

        # ----------------------------------------------------
        # Kartesische Koordinaten
        # ----------------------------------------------------

        x, y, z = earth_to_body

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
            z / distance_earth
        )

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

        distance_sun = np.linalg.norm(
            state_sun[:3]
        )

        velocity_sun = np.linalg.norm(
            state_sun[3:6]
        )

        # ----------------------------------------------------
        # Werte speichern
        # ----------------------------------------------------

        positions[body_name] = {
            "ra": ra,
            "dec": dec,
            "distance": distance_earth,
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
    f"Sonnensystem – Bezugspunkt Erdmittelpunkt – "
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

    if body_name == "Erde":

        # ----------------------------------------------------
        # Erde selbst hat vom Erdmittelpunkt aus keine
        # RA/Dec-Koordinate.
        #
        # Für die Darstellung wird deshalb ihre tatsächliche
        # Position auf der Umlaufbahn um die Sonne verwendet.
        # ----------------------------------------------------

        sun_state, _ = spice.spkezr(
            "SUN",
            start_et,
            "J2000",
            "NONE",
            "EARTH"
        )

        earth_to_sun = sun_state[:3]

        earth_position = -earth_to_sun

        distance = np.linalg.norm(
            earth_position
        )

        x, y, z = earth_position

        ra = np.arctan2(y, x)

        ra = (
            (ra + np.pi) % (2 * np.pi)
            - np.pi
        )

        dec = np.arcsin(
            z / distance
        )

    else:

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
# Infobox
# ============================================================

info_lines = []

for body_name in SOLSYS_DICT:

    data = positions[body_name]

    if np.isnan(data["ra"]):

        ra_text = "   ---"
        dec_text = "  ---"

    else:

        ra_deg = np.degrees(data["ra"])
        dec_deg = np.degrees(data["dec"])

        ra_text = f"{ra_deg:7.2f}°"
        dec_text = f"{dec_deg:6.2f}°"

        info_lines.append(
        f"{body_name.capitalize():8s} "
        f"{'Erde':4s} {data['distance_earth']/1e6:7.2f} Mio km   "
        f"{'Sonne':5s} {data['distance_sun']/1e6:7.2f} Mio km   "
        f"{'v':1s} {data['velocity_sun']:5.2f} km/s   "
        f"{'RA':2s} {ra_text}   "
        f"{'Dec':3s} {dec_text}"
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
    bbox_to_anchor=(1.25, 1.02)
)


# ============================================================
# Scaling / Layout
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
    # Planeten aktualisieren
    # --------------------------------------------------------

    for body_name in SOLSYS_DICT:

        if body_name == "Erde":

            sun_state, _ = spice.spkezr(
                "SUN",
                current_et,
                "J2000",
                "NONE",
                "EARTH"
            )

            earth_to_sun = sun_state[:3]

            earth_position = -earth_to_sun

            distance = np.linalg.norm(
                earth_position
            )

            x, y, z = earth_position

            ra = np.arctan2(y, x)

            ra = (
                (ra + np.pi) % (2 * np.pi)
                - np.pi
            )

            dec = np.arcsin(
                z / distance
            )

        else:

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
        f"Sonnensystem – Bezugspunkt Erdmittelpunkt – "
        f"{current_datetime.strftime('%Y-%m-%d %H:%M:%S UTC')}"
    )

    # --------------------------------------------------------
    # Infobox aktualisieren
    # --------------------------------------------------------

    info_lines = []

    for body_name in SOLSYS_DICT:

        data = positions[body_name]

        if np.isnan(data["ra"]):

            ra_text = f"{'---':>8s}"
            dec_text = f"{'---':>7s}"

        else:

            ra_deg = np.degrees(data["ra"])
            dec_deg = np.degrees(data["dec"])

            ra_text = f"{ra_deg:7.2f}°"
            dec_text = f"{dec_deg:6.2f}°"

        info_lines.append(
            f"{body_name.capitalize():8s} "
            f"{'Erde':4s} {data['distance_earth']/1e6:7.2f} Mio km   "
            f"{'Sonne':5s} {data['distance_sun']/1e6:7.2f} Mio km   "
            f"{'v':1s} {data['velocity_sun']:5.2f} km/s   "
            f"{'RA':2s} {ra_text}   "
            f"{'Dec':3s} {dec_text}"
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
