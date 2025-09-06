# Next / in progress:
## 0.5.6 Fitting TOI-4504

- welche der Auswertungen sind langsam? Diese seltener machen?

- lower und upper in mcmc result skalieren (wie die anderen Werte auch)

### TT-MCMC
- Residen fuer MaxLikelihood automatisch ausrechnen (Auch mean und median?)

- Eine Simulation laufen lassen. Deren entsprechende 12 c und 3 d Transits in 
tt_sim.csv
speichern und mcmc darauf laufen lassen.
- Wenn das konvergiert, sehe ich das als endgueltigen Beweis fuer Planet e!!!
- Mit 11 (ohne r, i), 13 (ohne r) und 15 Parametern testen. Evtl. mit 
  weiteren Parameterteilmengen.

- find_tts:
    - convert tts into a pandas Dataframe with columns eclipser, eclipsee, tt

- zunaechst nur moeglich nach TT ODER flux zu fitten
- spaeter gerne auch in Kombination

### Sonstige
- In mcmc Results: Max likelihood simflux vs measure flux plotten.  Oder 
direkt nur die Residuen plotten.
- RV-MCMC
    - zunaechst nur moeglich nach TT ODER flux ODER RV zu fitten
- Testen ob Performance besser, wenn bodies.find_tts einen Dataframe returniert
und CurveSimMCMC.residuals_tt_sum_squared entsprechend angepasst wird.
(Vermutlich nicht schneller, aber etwas eleganter?)

- Bei lokalen Minima:
    Separate MCMC-Laeufe mit stark eingeschraenktem lower-upper-Intervall?

- Nach jedem Chunk 5 Sekunden Zeit, aus einem Menue auszuwaehlen, 
  - sonst  laueft der naechste Chunk
  - Eine Option: Erstelle Video/Resultfile/Simflux fuer die aktuelle max 
    likelihood

- MCMC mit 20k Walker, Startwerte ueber sehr grossen Wertebereich fast 
gleichverteilt
- Andere Moves probieren (Parameter a von Stretching erhoehen brachte erstmal 
  nix)

- mcmc objekt abspeichern (pickle?) und spaeter (mit ggf. geandertem python 
code) fortsetzen

-Kompliziertere Fitting-Parameter Bereiche zulassen. Z.B. i in [86;88] oder 
[92;94]

### Housekeeping:
- Klasse in flux_data anlegen? (import * ist haesslich)

- Simulationscheck klappt fuer Drehimpuls, Energie und Center of Mass, aber 
nicht fuer Impuls.
    - Vielleicht Berechnungsfehler des Impulses?

- Lassen sich die Flux-Daten durch eine manuell von mir gemachte alternative 
Pixelmaske verbessern? Z.B. 4x4 statt 3x2?


## Programming Hinter MCMC-Integration zurueckgestellt:

### Result file:
- In save_mcmc_results() die "MCMC parameters" und "Program Parameters" 
sinnvoll vereinigen oder mit einer Hierarchie versehen.
- Irrelevante Attribute aus den Bodies entfernen, bevor Result-JSON 
  gespeichert wird
- **BUG!** Ich speichere fuer den Stern ueberfluessige Parameter mit falschen 
  Werten ab? Z.B. e und i.
- Orbitparameter regelmaessig (z.B. bei jedem Transit) im Resultfile speichern.


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
    - Extraplot im Video oder als PNG-file speichern: Lightcurve und die 
  Geschwindigkeit in Beobachterrichtung (y-Achse) für jeden Stern.
    - Plot radial velocity (movement of star to/from viewer). Provides a 
      check if the simulation is realistic.
    - Farbe des Sterns abhängig von RV varieren?

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
