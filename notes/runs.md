### Simulierte TT mit LMfit fitten
- Generiere mit den X024ML Parametern Simulierten Flux und TT 
- TOI4504c Inklination von 89.69 auf 89.79, damit Planet c nicht zu 
  frueh rausrutscht (also aufhoert Transits zu erzeugen)
- `TOI-4504_X024_maxl_dt20000.json` erzeugt
- Daraus mit results2list() und df2csv() in result_evaluation.py erzeugt:
  - `TT_TOI4504d.csv`
  - `TT_TOI4504c.csv`
- Daraus die TT rausgesucht, die ungefaehr den gemessenen TT (13c, 4d) 
  entsprechen und diese in `../data/simulations/tt_simX024.csv` gespeichert
- Configfile `../data/simulations/TOI-4504_simX024_01.ini` aus `TOI-4504_X024_maxl.ini` erzeugt
  - zunaechst die Body-Parameter nicht veraendert. D.h., LMfit 
  bekommt die optimalen Parameter-Werte bereits als Startwerte 
  serviert.
  - Nur 9 Parameter sind zu fitten:
    - d: m e P   ma
    - c: m e P O ma
  - Methoden powell, nelder, lbfgsb, cg, cobyla, bfgs, tnc, trust-constr, 
  ampgo:
    - starten und enden alle mit max_delta=0.0001, mean_delta=0.0000
  - Method differential_evolution divergiert trotz optimaler Startwerte
 
#### Konvergiert LMFit mit Startwerten in der Naehe (<1%) der richtigen Werte?
  - `TOI-4504_simX024_01.ini`
  - Methode nelder findet den ersten TT (direkt einen Tag nach StartDate) nicht
  - Methode powell findet den ersten TT nach ca. 60 Iterationen, aber kommt 
    trotzdem nicht nahe an die Parameter, mit denen die Daten erzeugt 
    wurden.
      - max_delta=1.4991, mean_delta=0.4837 [days] 

#### Konvergiert LMFit wenn ich mehr Parameter festhalte?
  - `TOI-4504_simX024_03.ini`
  - d.mass fitten, alle anderen Parameter auf 
    richtigen Wert festgesetzt
    - konvergiert _fuer jeden Startwert_ sehr schnell zu Residuen von cirka 
      null.
  - d.mass, d.e fitten: konvergiert zuverlaessig.
  - d.mass, d.e, d.P fitten: konvergiert ziemlich zuverlaessig, wenn der 
    Startwert fuer P weniger als 1 Tag vom richtigen Wert wegliegt.
  - d.mass, d.e, d.O, d.o, d.ma: konvergiert nicht
  - d.P, d.O, d.o, d.ma: konvergiert nicht
  - d.P, c.P: konvergiert nicht einmal mit brute force!

#### Konvergiert LMFit wenn ich 11 Parameter fitte mit sehr engen bounds?
- X024 ist ein MCMC fit mit 11 Parametern auf die 15 TT bis Sektor 94 (12c, 3d)
  - d: m e P O o ma
  - c: m e P _ o ma
  - MCMC ist ziemlich gut konvergiert, hat aber zu grosse Residuen
  - Ich fitte alle 11 Params mit LMfit (powell) auf die 17 TT inkl. Sektor 95 
    (13+4)
  - Dabei nutze ich die MCMC-Histogramme, um den Parametern sehr enge bounds 
    zu setzen
  - Konvergiert, aber zu grosse Residuen:
      - max_delta=1.4086   mean_delta=0.3836    [days] 
      - TOI4504d.mass:     1.27061
      - TOI4504d.e:        0.10229
      - TOI4504d.P:       41.16460
      - TOI4504d.Omega:   -0.14345
      - TOI4504d.omega:   88.10352
      - TOI4504d.ma:     328.52317
      - TOI4504c.mass:     3.63195
      - TOI4504c.e:        0.05419
      - TOI4504c.P:       83.28688
      - TOI4504c.omega:  295.95378
      - TOI4504c.ma:     151.22450
  - Moeglich Gruende fuer die zu grossen Residuen
      - LMfit konvergiert schlecht
      - Die Sektor 95-Daten passen schlecht zu den anderen Daten
      - Planet e
      - Fehler in CurveSimulator
  - Handlungsoptionen:
      - MCMC mit 11 Params auf den 13+4 TT laufen lassen. Danach evtl. 
        wieder LMfit mit MCMCs MaxL-Params als Startwerte nachlaufen lassen.
      - ???
 
- X026, X027 und X028 sind MCMC fits mit identischen Body- und Programmparametern
  - 11 Parametern auf die 17 TT bis Sektor 95 (13c, 4d) sind zu bestimmen
    - d: m e P O o ma
    - c: m e P _ o ma
  - **Residuen 15 std** nach ca. 4k steps incl. burn-in.
  - Die 3 Laeufe dienten dem Check, wie gut die Ergebnisse wiederholbar sind. 
    Nicht gut :(
  - Max Likelihood Parameter von X027 und X028 sind allerdings fast identisch.

- X029 ist wie X026 bis X028 aber mit c.P in 81.7-82.0 statt 82.0-83.7
  - Habe ich die ganze Zeit einen Fehler von Vitkova uebernommen und mit 
      falschen Startwerten/Bounds fuer P von c gearbeitet?
  - ausserdem bounds von omega und ma auf 0 bis 360 gesetzt
  - **Residuen 54 std** nach ca. 3k steps incl. burn-in. Enttaeuschend.
 
- X030
  - Erzeuge mit `TOI-4504_simX030_01.ini` `X030_sim_single_planet.json` 
    und daraus `tt_simX030.csv` mit 13 TT ungefaehr zeitlich so verteilt wie 
    die TESS Messungen von TOI4504c.
  - Veraendere Startwerte in ini-File leicht und starte LMfit
    - konvergiert schnell bei Residuen = 0
  - Bei stark veraenderten Startwerten stoppt LMfit mit powell oder nelder 
    schnell in einem lokalen Minimum.
  - Mit differential_evolution kommen die **Residuen** nach 11000 Iterationen 
    praktisch auf **Null**.
  - Dabei treffen erwartungsgemaess P genau, e groessenordnungsmaessig und die 
    anderen Parameter ueberhaupt nicht die Simulationsparameter.
  - MCMC verhaelt sich aehnlich: P wird super exakt bestimmt, e hat mean und 
    median in der Naehe, die anderen sind breit verteilt.
  - Deshalb jetzt X031!

- X031_c1d1_TT17
  - Nur c.P und d.P werden in diesem MCMC Run auf TT17 bestimmt!
  - Alle anderen Werte haben im wesentlichen die Vitkova-Werte.
  - **Residuen 311 std** nach 1250 steps incl. burn-in.
  - Acceptance 12%
 
- X032_c2d2_TT17
  - c.P, c.ma, d.P, d.ma werden in diesem MCMC Run auf TT17 bestimmt
  - Alle anderen Werte haben im wesentlichen die Vitkova-Werte.
  - **Residuen 169 std** nach 1750 steps incl. burn-in.
  - Acceptance 4%

- X033_c2d2e7_TT17
  - Wie X032 c.P, c.ma, d.P, d.ma aber zusaetzlich mit Planet e.meiPoom zu 
    fitten!
  - **Residuen 60 std** nach 5k steps incl. burn-in.
  - Acceptance 1.5%
 
- X034_c4d4e7_TT17
  - Wie X033 zusaetzlich c.m, c.e, d.m, d.e
  - **Residuen 45 std** nach 5k steps incl. burn-in.
  - Acceptance 2%

- X035_c7d7e7_TT17
  - Max Freiheit! Je 7 Params meiPOom fuer c, d und e.
  - **Residuen 30 std** nach 3250 steps incl. burn-in.
  - runs on Aspit, startet Sunday 19 Uhr

- X036_cd11P_TT17
  - wie X027, X028, aber
    - groesserer gueltiger Bereich fuer c.P
 
- X037_cd11P_TT17
  - wie X036, aber
    - noch groesserer gueltiger Bereich fuer c.P
    - -180 < ma < 180

#### Fits, bei denen c oder d oder e i>90 Grad hat

#### Vitkova-Fit (Nur 11 c transits bis Sektor 67)
...

#### Konvergiert LMFit fuer ein 1-Planeten-System?
...

#### Konvergiert LMFit fuer b (fest) und c (fitten), aber ohne d?
...


### MCMC fit von Flux inklusive Sektor 95
...

### MCMC fit von KEPLER-9
...

### LMfit fit von KEPLER-9
...

### Simulierten Flux mit LMfit fitten
- erfordert einige Programmier-Arbeit
 

