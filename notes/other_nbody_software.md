## Andere Software

### Exo Striker
- On Notebook2013, in Terminal, in directory ~/exostriker/exostriker, call 
python3.12 gui.py
- Unverstaendliche Fehlermeldungen
- letztes Update 2026
#### Fazit:
- Viele manuelle Eingriffe, bis es immerhin auf Linux laeuft
- Eingaben werden nur teilweise gespeicherrt (in Session)
- Kaum Dokumentation
- Ich krieg es nicht genutzt

### Juliet
- Doku: https://juliet.readthedocs.io/en/latest/
- Code: https://github.com/nespinoza/juliet
- Sprache: Python
- letztes Update 2026
#### Fazit:
- Kann kein n-body

### jnkepler
- Doku: https://jnkepler.readthedocs.io/en/stable/
- Code: https://github.com/kemasuda/jnkepler
- Sprache: Python
- letztes Update: 2026 
#### Fazit:
- Benutzt nur die abgewandelten Kepler-Parameter ecosw und esinw
- Keine Funktion fuer Single-Run dokumentiert? Beispielcode laesst direkt 
  einen Fit laufen.

### PlanetPlanet
- Code: https://github.com/rodluger/planetplanet
- Paper: https://iopscience.iop.org/article/10.3847/1538-4357/aa9c43/pdf
- Sprache: Python
- letztes Update 2017
#### Fazit:
- Mit viel Aufwand und AI-Hilfe uralten Code zum Laufen gebracht.
- Laut AI nutzt PlanetPlanet Rebound und will Jacobi Koordinaten. Ist also 
  eigentlich sehr geeignet fuer den Vergleich mit CurveSimulator.
- Die errechneten Transit Times sind aber hochgradig unglaubwuerdig. Der 
  erste TT con Planet c und d passt noch ungefaehr aber danach betraegt der 
  Abstand zwischen TT z.B. bei Planet d 43.5, 41.7, 43.1, 40.5, 40.5, 39.7, 39.5, 39.5, 39.2
- Ausserdem erzeugt Planet d angeblich Transits auch bei Start-Inklination 
  91.1.
- Gut moeglich, dass ich verkehrte Einheiten insbesondere der Abstaende oder 
  Massen benutzt habe.

### TTV Fast
- Doku+Code: https://github.com/simonrw/ttvfast-python 
- Sprache: Python (wrapper)
- letztes Update: 2018
#### Fazit:
- Fast TTV reproduziert die beobachteten TT mit Trifons Bet Fit Parametern sehr gut.
- Anpassung von Periode c alleine reicht nicht, um mit Ulis Best Fit die 
  beobachteten TT zu reproduzieren!
- Vielleicht in CurevSimulator Rebound durch ttvfast ersetzen?
- Nachteil: Nur 1-Stern-Systeme
- Nachteil: Nur TT, keine Info zu T1, T2, T3, T4
- Nachteil: Keine xyz-Koordinaten. Dadurch wohl keine Videos moeglich.
- Vorteil: Angeblich 5 bis 20 mal schneller als rebound. Vielleicht in MCMC 
  einsetzen?

#### exoplanet
- Doku: https://docs.exoplanet.codes/en/v0.1.4/
- Code: https://github.com/exoplanet-dev/exoplanet
- Sprache: Native Python Library
- Beherrscht sogar Dinge wie "Gaussian process models for stellar variability"
- Macht MCMC
- letztes Update: 2025 (last year) 

### PyDynamicaLC
- Code: https://github.com/YoffeG/PyDynamicaLC
- Sprache: Python
- letztes Update 2024
#### Fazit:
- ???

### Photodynamics.jl
- Paper: https://academic.oup.com/mnras/article/540/1/106/8125473?login=false
  - Table 1 listet diverse andere Softwarepakete auf
- Code: https://github.com/langfzac/Photodynamics.jl
- Doku: https://langfzac.github.io/Photodynamics.jl/dev/
- Sprache: Julia
- letztes Update 2025
#### Fazit:
- ???

### PhoDyMM
- Doku: 
- Code: https://github.com/dragozzine/PhoDyMM
- Paper1: https://ui.adsabs.harvard.edu/abs/2022BAAS...54e.360J/abstract
- Paper2: https://baas.aas.org/pub/2024n8i106p01/release/1
- Sprache: C (Python intern?)
- Betriebssystem: Linux
- letztes Update: 2025
#### Fazit:
- ?

### photodynam
- Paper: https://ui.adsabs.harvard.edu/abs/2017ascl.soft12013C/abstract
- Code: https://github.com/dfm/photodynam
- Doku: ???
- Sprache: C++
- letztes Update: 2013
#### Fazit:
- ?
