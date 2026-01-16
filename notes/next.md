# Next / in progress:
## 0.5.7 Automated Flux Data Download

### jetzt
- body.load() vereinfachen?
- complete make_max_likelihood_config_file() in mcmc
- Result json von F102 in bodies konvertieren und single run laufen lassen
  - flux Daten mit an die verschobenen Transits angepasstem tt-file neue 
    prozessieren
  - In den Einzel-Transits-Plots nach d-Transits suchen
  - flux_err_local statt flux_err benutzen
    - wie veraendert sich chi squared (plot)?
- warum sieht rv_o_vs_c.png inzwischen in der Mitte so komisch aus?
  - zumindest wenn aus mcmc heraus erzeugt, Beispiel F101 
  - wegen den zusaetzlichen Flux-Daten von theoretischen d-transits???

### aktuell
- Neue Kategorie Data in Configfile
 
- single_run macht eigentlich auch plots fuer die einzelnen 
  transits. Aber nicht beim mcmc single run! Vielleicht, weil ich tt-file 
  auskommentieren musste? Falls das der Grund war: moeglichkeit finden, dass 
  tt-file drin bleiben kann, obwohl mcmc nur mit flux daten arbeitet. 
  Differenzierteres Configfile erforderlich.


- single_run: Fallunterscheidungen sauberer machen
  - Je nach Szenario ist es unnoetig, dass die Simulation fuer 
    Standardabstaende dt berechnet wird (nur fuer Video sinnvoll?) 
  - Siehe # HACK: Hier wird body.positions nachtraeglich vergroessert, fuer 
    den Fall dass flux_file mehr Datenpunkte hat, als die Simulation

- neuer Parameter run_id
  - wird als unterverzeichnisname (statt laufender Nummer? oder neue 
    laufende Nummer hierrunter?) fuer die Ergebnisse benutzt
  - wird in (manchen) Dateinamen des Result-Outputs des Run benutzt
    - z.B. fuer gespeicherte Bodies Beispiel runid.TOI4504d.bdy


- Die einzelnen TT selber ordentlich fitten statt ExcelEyeballing
 
- WHFast convergence issue debuggen (warnung ist deaktiviert)

### Auswertungen 
NEXT! (vermutlich in cs_results.py reinschreiben)
#### Single run flux Auswertungen ausbauen und verbessern
- Fuer jeden Transit:
  - Zeitraum T1 bis T4 und +-x plotten (x = z.B. 2 Stunden).  
    - Einzelne Messungen (nur 120sek oder sogar alle Belichtungszeiten, jede 
      eine eigene Farbe)
    - binned flux (mehrere Linien mit unterschiedlicher Bin-Groesse)
    - computed flux
     
- tt mcmc erzeugt jetzt tt_delta.png UND tt_o_vs_c.png und die unterscheiden 
  sich auch noch geringfuegig!?!
- Ebenso Dopplung TOI-4504_T100.json und mcmc_results.json!!!

- Welche Single-Run-Auswertungen fehlen noch bei z.B. mcmc?

- measured_tt in Single-Run JSON-Result aufnehmen
   
- ChiSquared Verlauf Plot fuer mcmc?

- Sensitivitaetsanalyse: Plotten, wie veraendert sich Chi2, wenn man alle 
  gefundenen Parameter bis auf einen auf ihren MaxL-Werten festhaelt und 
  den letzten veraendert?

### Single Run 
- Hat single run jetzt alle plots/Auswertungen?

### MCMC
- nach einem Chunk auch die aktuelle Uhrzeit in die Konsole printen
- RV, TT, Flux und 
- Verschiedene Moves ausprobieren: 
  - https://emcee.readthedocs.io/en/stable/user/moves/

### Flux-MCMC
- mcmc flux wieder zum Laufen bringen
  - F056
  - laeuft jetzt, aber mean, maxL und median bleiben manchmal identisch und 
    konstant???
  - bei Gelegenheit debuggen
  
### LMfit
- Multi LMFit: normalverteilte Startwerte
- LMfit mit Fluxdaten laufen lassen? (erfordert neuen Code)

### TT-MCMC
- BUG: Nach Wiederaufnahme von MCMC mit gespeicherten Chains (.h5 file) starten 
  die Chains derzeit wieder bei den urspruenglichen Startwerten
- Nach jeder Verbesserung der MaxL, den Parametersatz abspeichern, damit ich 
  die Entwicklung der einzelnen MaxL-Parameter verfolgen kann

- find_tts:
    - convert tts into a pandas Dataframe with columns eclipser, eclipsee, tt

- Derzeit nur moeglich nach TT ODER flux zu fitten.
  - Kombination sinnvoll?

## Mit GUI (manuell) fitten/minimieren
- UPDATE PLOT ist noch buggy
- Vielleicht lieber selber eine Liste mit den 4 lines anlegen
  - schon bei einrichten des Plots
  - dann immer manuell die aelteste entfernen und die neueste hinzufuegen
- guifit.save_lmfit_results tut noch nix

### Sonstige
- Moeglichkeit entfernen, unterschiedliche dt konfigurieren zu koennen.
- In mcmc Results: Max likelihood simflux vs measure flux plotten.  Oder 
direkt nur die Residuen plotten.
- Testen ob Performance besser, wenn bodies.find_tts einen Dataframe returniert
und CurveSimMCMC.residuals_tt_sum_squared entsprechend angepasst wird.
(Vermutlich nicht schneller, aber etwas eleganter?)

- In einem File command.txt kann ich Anweisungen speichern, die nach jedem 
  Chunk oder sogar nach jeder Iteration ausgelesen werden. 
  - Eine Option: Erstelle Video/Resultfile/Simflux fuer die aktuelle max 
    likelihood
  - Erstelle jetzt saemtliche Result-Plots (statt nach jedem Chunk?)

- Andere emcee-Moves probieren (Parameter a von Stretching erhoehen brachte 
  erstmal 
  nix)

-Kompliziertere Fitting-Parameter Bereiche zulassen. 
  - Z.B. i in [86;88] oder [92;94]
  - Abhaengigkeit von anderen Fitting-Parametern

### RV-MCMC
- zunaechst nur moeglich nach TT ODER flux ODER RV zu fitten

### Housekeeping:
- checken ob die Bodynamen in eclipsers_names, eclipsees_names, primary 
  existieren

- Klasse in flux_data anlegen? (import * ist haesslich)

- Simulationscheck klappt fuer Drehimpuls, Energie und Center of Mass, aber 
nicht fuer Impuls.
    - Vielleicht Berechnungsfehler des Impulses?


## Programming Hinter MCMC-Integration zurueckgestellt:

### Result file:
- In save_mcmc_results() die "MCMC parameters" und "Program Parameters" 
sinnvoll vereinigen oder mit einer Hierarchie versehen.
- Irrelevante Attribute aus den Bodies entfernen, bevor Result-JSON 
  gespeichert wird
- **BUG!** Ich speichere fuer den Stern ueberfluessige Parameter mit falschen 
  Werten ab? Z.B. e und i.
- Orbitparameter regelmaessig (z.B. bei jedem Transit) im Resultfile speichern.

### Parameters/Config File
- Falls ich LMfit in CurveSimulator behalte: Den User warnen, wenn ein 
  illegaler Bodyname gewaehlt wurde (Kein Leerzeichen, kein Bindestrich, 
  keine fuehrende Ziffer)
- Extra Spalte bei Body-Params in Configfile mit dem Wert n oder u
  - n normal distribution (Spalte sigma ist std einer Gaussglocke)
  - u uniform distribution (Gleichverteilung von Startwert - sigma bis 
    Startwert + sigma)
 
### GUI:
- GUI fuer alles!!!
- Config File (optional) per GUI erstellen?
- Ermoeglichen, dass Body-Parameter von anderen Body-Parametern abhaengen?
    Bsp.: Radius planet c  = 1.2 * Radius planet d
- Add Neptun's mass and radius to astronomical constants in config file?
- Dem User selbstdefinierte Konstanten ermoeglichen, die dann im GUI benutzt 
  werden koennen?
    - Z.B. Radius Mars oder Luminosity Alpha Centauri
- Besser als eclipsers_names/eclipsees_names:
    - Jeder Body kriegt die Parameter relevant_eclipser  (yes/no) und 
  relevant_eclipsee (yes/no) zum Ankreuzen.

### Data Processing mit GUI
- Im GUI kann man auf Flux-Plots der einzelnen Sektoren mit der Maus T1, T2, TT,
  T3, T4 einzeichnen.
- Jeweils als Bereich (waagerechte Linie)
- Halbe Linienlaenge ist Standardabweichung
- Man zeichnet auch T0 und T5 ein
  - Die Intervalle T0:T1 und T4:T5 dienen zur Bestimmung des Median fuer 
    Normalisierung und Berechnung des flux_err (Standardabweichung der 
    Flux-Messungen)
- Man kann auch upper und lower limit einzeichnen, um outlier zu entfernen.

### Rebound:
- init_rebound() verallgemeinern. Auch Hierarchie ermoeglichen (=wer kreist 
um wen?)
- Systeme mit mehr als einem Stern korrekt hantieren
- Damit CurveSimulator auch fuer andere Systeme ausser TOI4504 funktioniert:
    - Aus dem User Input die Kepler-Parameter berechnen, die Rebound braucht, 
  und Rebound damit fuettern.
- Code entsorgen, der durch rebound abgedeckt wird.


### Configfile:
- Neu eingefuehrte Konstanten in Configfile statt hart im Code, z.B.
    - Mindestgenauigkeit bei binary search von Transittimes
    - Akzeptable Energieveraenderung
    - credible_mass
- Jeder Parameter soll seine eigene Skalierung haben koennen.
  - Also z.B. Sonnen-Masse anderen Faktor als Planetenmasse


### Result evaluation:
- komfortabel fuer User nutzbar machen.
    - z.B. Plots fuer Entwicklung von T14 oder b


### Physics:
- Meine Untersuchungen zum optimalen dt fuer ein konkretes System (Tabelle 
Rebound in TOI-4504.xlsx):
    - Kann ich das userfreundlich automatisieren?
    - Sonst zumindest halbautomatisch unterstuetzen und eine Anleitung fuer die Durchfuehrung schreiben.
- Multiple transits (und secondary eclipses?):
    - e.g. Triply Eclipsing Triples, >=3 bodies in line of sight
    - Handle them correctly
    - Example systems: https://arxiv.org/pdf/2204.02539.pdf
    - Siehe auch circle_eclipses.py?
- Indirect light / Planet Albedo
    - Bei Total System Luminosity beruecksichtigen -> kleiner dip in lightcurve bei secondary eclipse
    - Phasen beruecksichtigen: Halbplanet, Vollplanet, ...
    - Lohnt es sich die Verdunklung durch secondary eclipses und eclipses unter Planeten (fx planet c verdeckt planet d) zu modellieren?
      - Bei Albedo von 50% muesste eine secondary eclipse halb so tief sein 
        wie eine primary eclipse, korrekt?  
            - Nein, der Effekt ist viel kleiner und kann vorerst 
        vernachlaessigt werden. Kann ich irgendwann spaeter mal einbauen.
    - AI-Anfrage:   

            Wir beobachten die Helligkeit eines viele Lichtjahre entfernten Sternsystems mit einem Stern und einem Exoplaneten.
            Die Blickrichtung von der Erde zu dem Stern bildet mit der Verbindung von dem Stern zu dem Exoplaneten einen Winkel alpha.
            Wenn alpha = 0 Grad, befindet sich der Exoplanet von der Erde aus gesehen genau vor dem Stern. Bei alpha=180 Grad befindet
            sich der Exoplanet genau hinter dem Stern. Bei alpha=90 Grad befindet 
            sich der Exoplanet seitlich vom Stern und es wird von der Erde aus gesehen genau die Haelfte des Exoplaneten von dem Stern angeleuchtet. Sei max_light die maximale Lichtmenge
            (Luminosity), die der Exoplanet von dem Stern Richtung Erde reflektieren kann, wenn die komplette der Erde zugewandte
            Haelfte des Exoplaneten von dem Stern angeleuchtet wird.  Schreibe eine Pythonfunktion, die abhaengig von alpha errechnet,
            wie viel Prozent von max_light von dem Exoplaneten Richtung Erde reflektiert wird.
            Model variable stars?

- Parameter für regelmäßige Helligkeitsänderungen: variation = [(amplitude, 
periode, phase(shift?))]
  - Das deckt auch Sunspots ab. Die wären sonst zu beschwerlich.
  - Evtl. ermoeglichen, als Parameter den Namen einer Funktion anzugeben, die 
    die Helligkeit des Sterns abhaengig vom Zeitpunkt (in BJD?) berechnet.
- How to model sun spots?
    - Nora said (in episode 37, at 12:00): Spots move away from the center. Why?
    - Nora said (in episode 37, at 14:00): Peaks in frequency spectrum are integer multiples of peak with highest power. (harmonics) Why?
    - Is it sufficient to just model the lowest frequency? How big is the part I would miss?
- How to model pulsations and the like (sine curve)?
    - Nora said (in episode 37, at 15:30): Peaks in frequency spectrum are random.
- Should flares be modeled or are they too unpredictable?


### Housekeeping:
- Alle kritischen Textmeldungen mit colorama einfaerben (siehe color.py in 
Ordner debug)

- Warnungen unterdruecken??? Nervt jedenfalls waehrend mcmc
  - Vielleicht ist die Loesung auch, einen eigenen mcmc Fortschrittscounter zu 
    machen      
 

              WARNING in function find_tt: Rebound integration results are possibly not accurate enough.
              Try again with half the overall iteration time step parameter 'dt'.

- Im Config file steuern, ob Histogramm densities und bin_edges im 
MCMC-Result-JSON gespeichert werden sollen
    - Ist derzeit auskommentiert in mcmc_histograms()

- Package mal wieder hochladen und mit pip testen.
- Dependencies checken
    - Sind __init__.py und setup.py aktuell?
    - Dafuer auch benutzen: https://github.com/lichtgestalter/curvesimulator/network/dependencies


### Video:
    
- **BUG!!!!!!!** Seit ich keine Sampling rate mehr als Parameter habe laueft 
  irgendwas 
schief
    Beispiel: Teilweise werden weniger frames erzeugt als im Config File angegeben

- Orbits einzeichnen?
    - Kleiner Punkt (nur 1 Pixel?) fuer jeden Frame fuer jeden Koerper
    - Option draw_orbits (on/off) (pro Video oder pro Koerper?)
    - Option orbit_color

- Radial velocity:
    - Farbe des Sterns abhängig von RV varieren
      - Suche in cs_animation.py nach "Example code for changing circle color 
        during animation"
    - nur 3 bis 4 statt 5 bis 10 y-labels (ist oft so gequetscht)

- Option einbauen, die ein png Bild der Lichtkurve (und RV-Kurve?) erzeugt, 
statt ein Video (oder zusaetzlich)

- Realtime zwischen 2 frames ausrechnen und direkt im GUI anzeigen (muss 
nicht in jedem Intervall gleich sein)

### Dokumentation:

- Alles zu MCMC und rebound fehlt noch.

- Example videos
    - Videos erstellen und (auf youtube? eher vimeo wg. Bildqualitaet) 
  veröffentlichen.
    - Auf Video-Seite im Wiki (letzter Satz) ungefaehr sowas ergaenzen:
        - Or watch this [version with awful compression artefacts]
      (https://youtu.be/7diyCCIqiIc), courtesy of youtube.
  - Sternsysteme aus weiteren Papern visualisieren.
      Zunaechst nur Systeme mit nur einem Stern

- Dokumentation (Wiki) aktualisieren.
    - Wenn man mehrere Intervalle in starts/ends hat:
        - Anzahl in starts/ends/dt muss gleich sein (sonst wird Rest ignoriert)
        - Wenn lightcurve doof aussieht liegt das wahrscheinlich daran, dass 
          eines der Intervalle mitten in einem Transit startet oder endet.
        - Man kann starts/ends/dt komplett weglassen (Hinweis auf default Werte)

    - Results Output auf eigener Seite im Wiki dokumentieren.
    - Im Wiki die Videos verlinken.
    - Evtl. restliche Formulierungen von AI verbessern lassen.

              Take a look at this text on a page of a wiki about a python library. The text is written in mark-up language.
              Improve the language, for example by removing spelling errors and making the formulations more elegant.
              Format your improvements in the same mark-up language: ""

- Um Hinweis auf CurveSimulator(-Paper?) bei Veroeffentlichungen bitten
    - Im GitHubWiki
    - Print in Konsole
    - Im JSON File
    - Copyright notice ins Video?

## Ferne Zukunft / Neues Projekt
### Einen manuellen/visuelle Fitter bauen
- Man waehlt 2 der Fittingparameter aus.
- Bekommt ein Diagramm der bisher ermittelten Residuen angezeigt
  - farbcodiert Z.B. je gruener, desto kleineres Residuum.
  - besser noch: Graustufen und Top n% (user definierbar) in rot
  - man kann waehlen ob die Farbcodierung global ist oder sich beim 
    Reinzoomen anpasst
  - vielleicht global Grautoene und lokal Blautoene nehmen
- Man kann in das Diagramm entweder eine Linie oder ein Gitter zeichnen.
- Daraufhin werden die Residuen fuer alle Kombinationen der beiden 
  Parameter entlang der eingezeichneten Linien berechnet und geplottet.
- Man kann den Wertebereich des Diagramms verschieben und hineinzoomen. 
- Wenn der Fitter nix zu tun hat, geht er zu lokalen Minima uns sucht drumherum.
- Man kann auch waehrend des Fits Startwerte definieren von denen aus mit 
  einer zu waehlenden LMfit-Methode gesucht wird. Auch die LMfit runs gehen 
  in die Datenbank aller Residuen ein.