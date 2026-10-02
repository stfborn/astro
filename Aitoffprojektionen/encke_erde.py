import re
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from datetime import datetime, timezone, timedelta
from pathlib import Path
import spiceypy as spice


# ============================================================
# Dateien
# Encke-Kernel mit encke_kernel_erstellen.py erzeugt
# (Horizons-Vektoren, Type 8, 2020–2036)
# ============================================================

BASE = Path(__file__).resolve().parent

LSK = BASE / "naif0012.tls"
SPK = BASE / "de440.bsp"
ENCK = BASE / "encke_2020_2036.bsp"

# NAIF-ID 2P/Encke (als Text, so will es spkezr)
ENCKE = "1000002"


# ============================================================
# Konstanten (nur für den kernel-losen Rückfallweg)
# ============================================================

AU_KM = 149597870.7
GM_SONNE = 1.32712440018e11      # km^3/s^2
K_GAUSS = 0.01720209895          # Gaußsche Konstante, Grad/Tag
EPSILON = np.radians(23.4392911) # Schiefe der Ekliptik (J2000)


# ============================================================
# 2P/Encke – Bahnelemente als Rückfallwert
# (nur wenn der Kernel fehlt)
# Quelle: JPL Horizons, 01.10.2026
# ============================================================

ENCKE_FALLBACK = {
    "epoche_jd": 2461315.214735802,
    "a_au": 2.217731397662715,
    "e": 0.8473172110473277,
    "i_deg": 11.34776328477018,
    "Om_deg": 334.0189642734944,
    "w_deg": 187.2876797339125,
    "M_deg": 320.7524895514863,
    "quelle": "eingebettet (Horizons-Epoche 01.10.2026)",
}


# ============================================================
# Bahnelemente laden (nur ohne Kernel nötig)
# ============================================================

def encke_elemente_laden():

    try:
        from astroquery.jplhorizons import Horizons
        from astropy.time import Time
    except ImportError:
        print("astroquery fehlt – nutze eingebettete Elemente.")
        return dict(ENCKE_FALLBACK)

    try:
        jd = float(Time.now().jd)

        try:
            tab = Horizons(
                id="2P",
                location="500@10",
                epochs=[jd],
                id_type="smallbody",
            ).elements()
        except ValueError as err:
            # Mehrere Erscheinungen -> höchste
            # Datensatznummer ist die aktuelle
            records = sorted(
                {int(x) for x in re.findall(r"90\d{6}", str(err))}
            )
            if not records:
                raise
            tab = Horizons(
                id=str(records[-1]),
                location="500@10",
                epochs=[jd],
            ).elements()

        reihe = tab[0]

        print("Encke-Elemente live von Horizons.")

        return {
            "epoche_jd": float(tab["datetime_jd"][0]),
            "a_au": float(reihe["a"]),
            "e": float(reihe["e"]),
            "i_deg": float(reihe["incl"]),
            "Om_deg": float(reihe["Omega"]),
            "w_deg": float(reihe["w"]),
            "M_deg": float(reihe["M"]),
            "quelle": "JPL Horizons live",
        }

    except Exception as err:
        print(f"Horizons-Abruf gescheitert ({err})")
        print("Nutze eingebettete Elemente.")
        return dict(ENCKE_FALLBACK)


# ============================================================
# SPICE Kernel laden
# ============================================================

print("LSK:", LSK)
print("LSK vorhanden:", LSK.exists())

print("SPK:", SPK)
print("SPK vorhanden:", SPK.exists())

print("ENCK:", ENCK)
print("ENCK vorhanden:", ENCK.exists())

spice.furnsh(str(LSK))
spice.furnsh(str(SPK))

KERNEL_VORHANDEN = ENCK.exists()

if KERNEL_VORHANDEN:
    spice.furnsh(str(ENCK))
    print("Encke aus Kernel – läuft offline.")
else:
    print("Encke-Kernel fehlt – Rückfall auf Kepler-Ausbreitung.")


# ============================================================
# Kontrolle
# ============================================================

print()
print("Geladene Kernel:", spice.ktotal("ALL"))

for i in range(spice.ktotal("ALL")):
    print(spice.kdata(i, "ALL"))


# ============================================================
# Startzeit + Elemente (nur ohne Kernel)
# ============================================================

start_datetime = datetime.now(timezone.utc)

start_et = spice.str2et(
    start_datetime.strftime("%Y-%m-%d %H:%M:%S UTC")
)

start_monotonic = time.monotonic()

if KERNEL_VORHANDEN:
    ELEMENTE = {"quelle": f"Kernel {ENCK.name} (offline)"}
else:
    ELEMENTE = encke_elemente_laden()

print()
print("Encke-Quelle:", ELEMENTE["quelle"])


# ============================================================
# Körper (SSB-Sicht wie in solar_system_barycenter.py)
# ============================================================

KOERPER = {
    "Encke": ENCKE,
    "Erde": "EARTH",
}

KOERPER_FARBE = {
    "Encke": "tab:cyan",
    "Erde": "tab:blue",
}


# ============================================================
# Kepler-Ausbreitung (nur ohne Kernel)
# ============================================================

def loese_kepler(M_rad, e):

    E = M_rad if e < 0.8 else np.pi

    for _ in range(60):
        delta = (
            (E - e * np.sin(E) - M_rad)
            / (1.0 - e * np.cos(E))
        )
        E -= delta
        if abs(delta) < 1e-13:
            break

    return E


def encke_position(elemente, jd):
    """
    Heliozentrische Encke-Position im J2000-Äquatorsystem.

    Gibt (Positionsvektor [km], Sonnenabstand [AE],
    Bahngeschwindigkeit [km/s]) zurück.
    """

    a = elemente["a_au"]
    e = elemente["e"]
    i = np.radians(elemente["i_deg"])
    Om = np.radians(elemente["Om_deg"])
    w = np.radians(elemente["w_deg"])

    # --------------------------------------------------------
    # Mittlere Anomalie zur Zielzeit
    # --------------------------------------------------------

    n = K_GAUSS / a ** 1.5

    M = np.radians(
        (elemente["M_deg"] + n * (jd - elemente["epoche_jd"])) % 360.0
    )

    # --------------------------------------------------------
    # Exzentrische -> wahre Anomalie
    # --------------------------------------------------------

    E = loese_kepler(M, e)

    nu = 2.0 * np.arctan2(
        np.sqrt(1.0 + e) * np.sin(E / 2.0),
        np.sqrt(1.0 - e) * np.cos(E / 2.0)
    )

    r_au = a * (1.0 - e * np.cos(E))

    # --------------------------------------------------------
    # Perifokal -> Ekliptik (J2000)
    # --------------------------------------------------------

    xw = r_au * np.cos(nu)
    yw = r_au * np.sin(nu)

    cosO, sinO = np.cos(Om), np.sin(Om)
    cosi, sini = np.cos(i), np.sin(i)
    cosw, sinw = np.cos(w), np.sin(w)

    x_ekl = (
        (cosO * cosw - sinO * sinw * cosi) * xw
        + (-cosO * sinw - sinO * cosw * cosi) * yw
    )
    y_ekl = (
        (sinO * cosw + cosO * sinw * cosi) * xw
        + (-sinO * sinw + cosO * cosw * cosi) * yw
    )
    z_ekl = (sinw * sini) * xw + (cosw * sini) * yw

    # --------------------------------------------------------
    # Ekliptik -> Äquatorial (J2000)
    # --------------------------------------------------------

    x = x_ekl
    y = np.cos(EPSILON) * y_ekl - np.sin(EPSILON) * z_ekl
    z = np.sin(EPSILON) * y_ekl + np.cos(EPSILON) * z_ekl

    position_km = np.array([x, y, z]) * AU_KM

    # --------------------------------------------------------
    # Bahngeschwindigkeit (Vis-Viva)
    # --------------------------------------------------------

    v = np.sqrt(
        GM_SONNE * (2.0 / (r_au * AU_KM) - 1.0 / (a * AU_KM))
    )

    return position_km, r_au, v


# ============================================================
# Positionen bestimmen
# ============================================================

def bahn_daten(name, et):
    """
    Gibt (Position SSB [km], Position Sonne-relativ [km],
    Position Erde-relativ [km], Geschwindigkeit [km/s]) zurück.
    """

    if name == "Encke" and not KERNEL_VORHANDEN:

        jd = 2451545.0 + et / 86400.0

        helio, _, v = encke_position(ELEMENTE, jd)

        sonne_ssb, _ = spice.spkezr(
            "SUN", et, "J2000", "NONE", "SSB"
        )

        erde_ssb, _ = spice.spkezr(
            "EARTH", et, "J2000", "NONE", "SSB"
        )

        erde_helio = erde_ssb[:3] - sonne_ssb[:3]

        return (
            sonne_ssb[:3] + helio,
            helio,
            helio - erde_helio,
            v,
        )

    target = KOERPER[name]

    zustand_ssb, _ = spice.spkezr(
        target,
        et,
        "J2000",
        "LT+S",
        "SSB"
    )

    zustand_sonne, _ = spice.spkezr(
        target,
        et,
        "J2000",
        "LT+S",
        "SUN"
    )

    zustand_erde, _ = spice.spkezr(
        target,
        et,
        "J2000",
        "LT+S",
        "EARTH"
    )

    return (
        zustand_ssb[:3],
        zustand_sonne[:3],
        zustand_erde[:3],
        float(np.linalg.norm(zustand_sonne[3:6])),
    )


def calculate_positions(et):

    positions = {}

    for name in KOERPER:

        pos_ssb, pos_sonne, pos_erde, v = bahn_daten(name, et)

        x, y, z = pos_ssb

        distance_ssb = float(np.linalg.norm(pos_ssb))
        distance_sun = float(np.linalg.norm(pos_sonne))
        distance_earth = float(np.linalg.norm(pos_erde))

        # ----------------------------------------------------
        # Rektaszension
        # --------------------------------------------------------

        ra = np.arctan2(y, x)

        ra = (ra + np.pi) % (2 * np.pi) - np.pi

        # ----------------------------------------------------
        # Deklination
        # --------------------------------------------------------

        dec = np.arcsin(z / distance_ssb)

        positions[name] = {
            "ra": ra,
            "dec": dec,
            "distance": distance_ssb,
            "distance_sun": distance_sun,
            "distance_earth": distance_earth,
            "velocity_sun": v,
        }

    return positions


def abstand_zeile(positions):
    d = positions["Encke"]["distance_earth"]
    return (
        f"2P/Encke – Erde: {d / 1e6:8.2f} Mio km  "
        f"({d / AU_KM:6.4f} AE)"
    )


def info_zeilen(positions):
    zeilen = []
    for name in KOERPER:
        daten = positions[name]
        zeilen.append(
            f"{name:6s} "
            f"SSB {daten['distance'] / 1e6:7.2f} Mio km   "
            f"Sonne {daten['distance_sun'] / 1e6:7.2f} Mio km   "
            f"Erde {daten['distance_earth'] / 1e6:7.2f} Mio km   "
            f"v {daten['velocity_sun']:5.2f} km/s   "
            f"RA {np.degrees(daten['ra']):7.2f}°   "
            f"Dec {np.degrees(daten['dec']):6.2f}°"
        )
    return zeilen


# ============================================================
# Erste Position
# ============================================================

positions = calculate_positions(start_et)

print()
print(abstand_zeile(positions))
print()


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
    f"2P/Encke – Erde (SSB) – "
    f"{start_datetime.strftime('%Y-%m-%d %H:%M:%S UTC')}",
    fontsize=12,
    fontweight="bold",
    x=0.50,
    y=1.06
)


# ============================================================
# Körper zeichnen
# ============================================================

punkte = {}

for name in KOERPER:

    punkt, = ax.plot(
        [positions[name]["ra"]],
        [positions[name]["dec"]],
        color=KOERPER_FARBE[name],
        marker="o",
        linestyle="None",
        markersize=10,
        label=name
    )

    punkte[name] = punkt


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
# Ekliptik
# ============================================================

ra_ecliptic = np.linspace(
    -np.pi,
    np.pi,
    1000
)

dec_ecliptic = np.arcsin(
    np.sin(EPSILON)
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

plt.legend(
    loc="upper right",
    bbox_to_anchor=(1.15, 1.02)
)


# ============================================================
# Infobox: Details
# ============================================================

info_text = fig.text(
    0.5,
    0.14,
    "\n".join(info_zeilen(positions)),
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
# Scaling / Layout
# ============================================================

plt.subplots_adjust(
    bottom=0.32,
    top=0.92,
    right=0.88
)


# ============================================================
# Animation (alle 2 Sekunden)
# ============================================================

def update(frame):

    # --------------------------------------------------------
    # Tatsächlich vergangene Zeit
    # --------------------------------------------------------
    # 1 reale Sekunde = 1 Sekunde Simulation

    elapsed_seconds = time.monotonic() - start_monotonic

    current_et = start_et + elapsed_seconds

    current_datetime = (
        start_datetime
        + timedelta(seconds=elapsed_seconds)
    )

    # --------------------------------------------------------
    # Neue Positionen + Entfernung
    # --------------------------------------------------------

    positions = calculate_positions(current_et)

    # --------------------------------------------------------
    # Punkte aktualisieren
    # --------------------------------------------------------

    for name in KOERPER:
        punkte[name].set_data(
            [positions[name]["ra"]],
            [positions[name]["dec"]]
        )

    # --------------------------------------------------------
    # Titel aktualisieren
    # --------------------------------------------------------

    title.set_text(
        f"2P/Encke – Erde (SSB) – "
        f"{current_datetime.strftime('%Y-%m-%d %H:%M:%S UTC')}"
    )

    # --------------------------------------------------------
    # Infobox aktualisieren
    # --------------------------------------------------------

    info_text.set_text(
        "\n".join(info_zeilen(positions))
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
