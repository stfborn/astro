# Planeten über Mainz

Vier eigenständige Skripte (Skyfield + JPL DE421, `de421.bsp` liegt im
selben Ordner). Standort-Standard: Mainz, 49.9858° N, 8.2791° E.
Alle Zeiten in Europe/Berlin (MESZ/MEZ).

## mars_mainz.py, venus_mainz.py, jupiter_mainz.py

Berechnen pro Tag: Aufgang / Kulmination (höchster Stand mit Höhe und
Azimut) / Untergang, dazu einen Tagesverlauf (Höhe, Azimut,
Himmelsrichtung, RA/Dec, Entfernung), einen „Jetzt“-Block und enge
Mond-Konjunktionen unter 5° — letztere in **Rot**.

```sh
python3 mars_mainz.py --help
python3 mars_mainz.py
python3 mars_mainz.py --datum 2026-09-28 --tage 7
python3 mars_mainz.py --verlauf 2026-09-28 --schritt 60
python3 jupiter_mainz.py --tage 14 --kein-verlauf
```

Wichtigste Optionen: `--datum JJJJ-MM-TT` (Standard: heute),
`--tage N` (Standard: 7), `--schritt MIN` (Standard: 60),
`--verlauf JJJJ-MM-TT`, `--kein-verlauf`, `--lat/--lon/--elev`.

## mond_planet_linien.py

Linien-Diagramm des Winkelabstands Mond–Mars, Mond–Venus und
Mond–Jupiter über einen Zeitraum. Wo eine Linie ein Minimum hat und
die Schwellenlinie unterschreitet, steht eine enge Begegnung an
(Markierung mit Datum). Die Ereignisliste erscheint zusätzlich im
Terminal als PNG-Datei.

```sh
python3 mond_planet_linien.py --help
python3 mond_planet_linien.py --datum 2026-09-28 --tage 30
python3 mond_planet_linien.py --tage 30 --schwelle 3 --output begegnungen.png
```

Wichtigste Optionen: `--tage N` (Standard: 30), `--schritt STD`
(Raster in Stunden, Standard: 3), `--schwelle GRAD` (Standard: 5),
`--output DATEI.png`.

Hinweis: Echte „Oppositionen“ (Sonne–Erde–Planet = 180°) gibt es nur
bei Mars und Jupiter — hier geht es um enge Mond–Planet-Begegnungen
(Konjunktionen) am Himmel.
