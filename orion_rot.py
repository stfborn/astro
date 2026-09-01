import numpy as np
import matplotlib.pyplot as plt

import os
import sys

from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider, Button

from datetime import datetime, timedelta, timezone

from skyfield.api import load, Star, wgs84


# ============================================================
# EINSTELLUNGEN
# ============================================================

LATITUDE = 49.9858
LONGITUDE = 8.2791
ELEVATION = 90

START = datetime(2026, 8, 15, 0, 0)
END = datetime(2026, 12, 31, 23, 0)

STEP = timedelta(minutes=15)


# ============================================================
# SOMMER-/WINTERZEIT
# ============================================================

def mainz_to_utc(local_time):

    if local_time < datetime(2026, 10, 25):

        offset = timedelta(hours=2)

    else:

        offset = timedelta(hours=1)

    return (local_time - offset).replace(
        tzinfo=timezone.utc
    )


# ============================================================
# SKYFIELD
# ============================================================

print()
print("Skyfield wird gestartet...")
print()

ts = load.timescale()

print("Lade JPL-Ephemeride DE421...")

if getattr(sys, "frozen", False):
    base_path = sys._MEIPASS
else:
    base_path = os.path.dirname(
        os.path.abspath(__file__)
    )

de421_path = os.path.join(
    base_path,
    "de421.bsp"
)

if not os.path.exists(de421_path):
    raise FileNotFoundError(
        f"DE421 nicht gefunden: {de421_path}"
    )

eph = load(de421_path)

earth = eph["earth"]

mainz = earth + wgs84.latlon(
    LATITUDE,
    LONGITUDE,
    elevation_m=ELEVATION
)

print("DE421 geladen.")
print()


# ============================================================
# ORION
# ============================================================

orion_stars = {

    "Betelgeuse":
        Star(
            ra_hours=5 + 55/60 + 10.3/3600,
            dec_degrees=7 + 24/60 + 25/3600
        ),

    "Bellatrix":
        Star(
            ra_hours=5 + 25/60 + 7.9/3600,
            dec_degrees=6 + 20/60 + 59/3600
        ),

    "Mintaka":
        Star(
            ra_hours=5 + 32/60 + 0.4/3600,
            dec_degrees=-(0 + 17/60 + 57/3600)
        ),

    "Alnilam":
        Star(
            ra_hours=5 + 36/60 + 14.5/3600,
            dec_degrees=-(1 + 12/60 + 7/3600)
        ),

    "Alnitak":
        Star(
            ra_hours=5 + 40/60 + 45.5/3600,
            dec_degrees=-(1 + 56/60 + 34/3600)
        ),

    "Saiph":
        Star(
            ra_hours=5 + 47/60 + 45.4/3600,
            dec_degrees=-(9 + 40/60 + 10/3600)
        ),

    "Rigel":
        Star(
            ra_hours=5 + 14/60 + 32.3/3600,
            dec_degrees=-(8 + 12/60 + 6/3600)
        ),

    "Meissa":
        Star(
            ra_hours=5 + 35/60 + 14.5/3600,
            dec_degrees=9 + 56/60 + 3/3600
        ),
}


# ============================================================
# ORION-LINIEN
# ============================================================

connections = [
    ("Betelgeuse", "Bellatrix"),
    ("Betelgeuse", "Meissa"),
    ("Meissa", "Bellatrix"),

    ("Bellatrix", "Mintaka"),
    ("Mintaka", "Alnilam"),
    ("Alnilam", "Alnitak"),

    ("Betelgeuse", "Alnitak"),
    ("Alnitak", "Saiph"),

    ("Mintaka", "Rigel"),
]

# ============================================================
# ZEITPUNKTE
# ============================================================

local_dates = []

current = START

while current <= END:

    local_dates.append(current)

    current += STEP


# ============================================================
# UTC
# ============================================================

utc_dates = [
    mainz_to_utc(t)
    for t in local_dates
]


# ============================================================
# SKYFIELD ZEITEN
# ============================================================

sky_times = ts.from_datetimes(
    utc_dates
)


# ============================================================
# POSITIONEN BERECHNEN
# ============================================================

positions = {}


for name, star in orion_stars.items():

    print(
        f"Berechne {name}..."
    )

    astrometric = mainz.at(
        sky_times
    ).observe(star)

    apparent = astrometric.apparent()

    altitude, azimuth, distance = apparent.altaz()

    positions[name] = {

        "altitude": altitude.degrees,

        "azimuth": azimuth.degrees
    }


print()
print("Berechnung abgeschlossen.")
print()


# ============================================================
# FIGURE
# ============================================================

fig = plt.figure(
    figsize=(11, 9)
)

# Sehr dunkler Hintergrund
fig.patch.set_facecolor("#050B18")

ax = fig.add_subplot(
    111
)

ax.set_facecolor("#050B18")


# ============================================================
# STERNE
# ============================================================

star_points = {}

for name in orion_stars:

    star_points[name], = ax.plot(
        [],
        [],
        "o",
        markersize=9,
	color='red'
    )


# ============================================================
# KONSTELLATIONSLINIEN
# ============================================================

connection_lines = []

for star1, star2 in connections:

    line, = ax.plot(
        [],
        [],
        "-",
        linewidth=1.2,
        color="red"
    )

    connection_lines.append(
        (star1, star2, line)
    )


# ============================================================
# STERNNAMEN
# ============================================================

labels = {}

for name in orion_stars:

    labels[name] = ax.text(
        0,
        0,
        name,
        fontsize=9,
        color="red"
    )


# ============================================================
# DATUM
# ============================================================

date_text = ax.text(
    0.5,
    1.04,
    "",
    transform=ax.transAxes,
    ha="center",
    fontsize=16,
    color="white"
)


# ============================================================
# INFO
# ============================================================

info_text = ax.text(
    0.5,
    1.005,
    "Mainz 55130 – Azimut / Höhe",
    transform=ax.transAxes,
    ha="center",
    fontsize=11,
    color="white"
)


# ============================================================
# ACHSEN
# ============================================================

ax.set_xlabel(
    "Azimut (Grad)",
    color="white"
)

ax.set_ylabel(
    "Höhe über Horizont (Grad)",
    color="white"
)


# Achsen und Beschriftungen hell
ax.tick_params(
    axis="both",
    colors="white"
)


# Rahmen
for spine in ax.spines.values():

    spine.set_color("#444444")


ax.set_xlim(
    0,
    360
)

ax.set_ylim(
    0,
    90
)


ax.set_xticks([
    0,
    45,
    90,
    135,
    180,
    225,
    270,
    315,
    360
])


ax.set_xticklabels([
    "N",
    "NO",
    "O",
    "SO",
    "S",
    "SW",
    "W",
    "NW",
    "N"
])


ax.grid(
    True,
    color="#30343b",
    alpha=0.5
)


# ============================================================
# ZOOM-SLIDER
# ============================================================

zoom_ax = plt.axes([
    0.18,
    0.085,
    0.64,
    0.025
])

zoom_slider = Slider(
    zoom_ax,
    "Zoom",
    0,
    1,
    valinit=0.92,
    valstep=0.01
)
zoom_slider.label.set_color("white")
zoom_slider.valtext.set_color("white")

# ============================================================
# ZEIT-SLIDER
# ============================================================

time_ax = plt.axes([
    0.18,
    0.045,
    0.64,
    0.025
])

initial_frame = 22

time_slider = Slider(
    time_ax,
    "Zeit",
    0,
    len(local_dates) - 1,
    valinit=initial_frame,
    valstep=1
)
time_slider.label.set_color("white")
time_slider.valtext.set_color("white")

# ============================================================
# BUTTONS
# ============================================================

back_ax = plt.axes([
    0.15,
    0.005,
    0.12,
    0.03
])

play_ax = plt.axes([
    0.44,
    0.005,
    0.12,
    0.03
])

forward_ax = plt.axes([
    0.73,
    0.005,
    0.12,
    0.03
])


back_button = Button(
    back_ax,
    "◀ Zurück"
)

play_button = Button(
    play_ax,
    "▶ Play"
)

forward_button = Button(
    forward_ax,
    "Weiter ▶"
)


# ============================================================
# STATUS
# ============================================================

current_frame = initial_frame

playing = False


# ============================================================
# ZOOM-BEREICH
# ============================================================

def update_zoom():

    zoom = zoom_slider.val

    # Gesamthimmel
    if zoom == 0:

        ax.set_xlim(
            0,
            360
        )

        ax.set_ylim(
            0,
            90
        )

        return


    frame = current_frame

    azimuths = []
    altitudes = []


    for name in orion_stars:

        altitude = positions[name][
            "altitude"
        ][frame]

        azimuth = positions[name][
            "azimuth"
        ][frame]


        if altitude >= 0:

            azimuths.append(
                azimuth
            )

            altitudes.append(
                altitude
            )


    if not azimuths:

        return


    center_az = np.mean(
        azimuths
    )

    center_alt = np.mean(
        altitudes
    )


    # Zoom
    width = 360 - zoom * 340

    height = 90 - zoom * 70


    x_min = center_az - width / 2

    x_max = center_az + width / 2

    y_min = max(
        0,
        center_alt - height / 2
    )

    y_max = min(
        90,
        center_alt + height / 2
    )


    if x_min < 0:

        x_min += 360
        x_max += 360


    if x_max > 360:

        x_min -= 360
        x_max -= 360


    ax.set_xlim(
        x_min,
        x_max
    )

    ax.set_ylim(
        y_min,
        y_max
    )


# ============================================================
# DARSTELLUNG AKTUALISIEREN
# ============================================================

def update(frame):

    global current_frame

    current_frame = int(frame)

    local_time = local_dates[
        current_frame
    ]


    # ========================================================
    # STERNE
    # ========================================================

    for name in orion_stars:

        altitude = positions[name][
            "altitude"
        ][current_frame]

        azimuth = positions[name][
            "azimuth"
        ][current_frame]


        if altitude < 0:

            star_points[name].set_data(
                [],
                []
            )

            labels[name].set_visible(
                False
            )

        else:

            star_points[name].set_data(
                [azimuth],
                [altitude]
            )

            labels[name].set_position(
                (
                    azimuth,
                    altitude + 1.5
                )
            )

            labels[name].set_visible(
                True
            )


    # ========================================================
    # ORION-LINIEN
    # ========================================================

    for star1, star2, line in connection_lines:

        alt1 = positions[star1][
            "altitude"
        ][current_frame]

        az1 = positions[star1][
            "azimuth"
        ][current_frame]

        alt2 = positions[star2][
            "altitude"
        ][current_frame]

        az2 = positions[star2][
            "azimuth"
        ][current_frame]


        if alt1 >= 0 and alt2 >= 0:

            line.set_data(
                [az1, az2],
                [alt1, alt2]
            )

            line.set_visible(
                True
            )

        else:

            line.set_data(
                [],
                []
            )

            line.set_visible(
                False
            )


    # ========================================================
    # DATUM
    # ========================================================

    date_text.set_text(
        local_time.strftime(
            "%d.%m.%Y – %H:%M Uhr"
        )
    )


    # ========================================================
    # ZOOM
    # ========================================================

    update_zoom()


    return (
        list(star_points.values())
        +
        [
            x[2]
            for x in connection_lines
        ]
        +
        list(labels.values())
        +
        [date_text]
    )


# ============================================================
# ZEIT-SLIDER
# ============================================================

def time_changed(value):

    update(
        int(value)
    )

    fig.canvas.draw_idle()


time_slider.on_changed(
    time_changed
)


# ============================================================
# ZOOM-SLIDER
# ============================================================

def zoom_changed(value):

    update_zoom()

    fig.canvas.draw_idle()


zoom_slider.on_changed(
    zoom_changed
)


# ============================================================
# ZURÜCK
# ============================================================

def go_back(event):

    global current_frame

    if current_frame > 0:

        current_frame -= 1

        time_slider.set_val(
            current_frame
        )


back_button.on_clicked(
    go_back
)


# ============================================================
# WEITER
# ============================================================

def go_forward(event):

    global current_frame

    if current_frame < len(local_dates) - 1:

        current_frame += 1

        time_slider.set_val(
            current_frame
        )


forward_button.on_clicked(
    go_forward
)


# ============================================================
# PLAY / PAUSE
# ============================================================

def toggle_play(event):

    global playing

    playing = not playing


    if playing:

        play_button.label.set_text(
            "⏸ Pause"
        )

        animation.event_source.start()

    else:

        play_button.label.set_text(
            "▶ Play"
        )

        animation.event_source.stop()


play_button.on_clicked(
    toggle_play
)


# ============================================================
# ANIMATION
# ============================================================

def animation_update(frame):

    global current_frame

    if not playing:

        return []


    current_frame += 1


    if current_frame >= len(local_dates):

        current_frame = 0


    update(
        current_frame
    )


    time_slider.set_val(
        current_frame
    )


    return []


animation = FuncAnimation(
    fig,
    animation_update,
    interval=80,
    blit=False,
    repeat=True
)


# ============================================================
# START
# ============================================================

update(
    initial_frame
)


plt.title(
    "Orion über Mainz – Zoomdarstellung",
    fontsize=18,
    color="white"
)


plt.subplots_adjust(
    bottom=0.16,
    top=0.88
)


plt.show()
