✅ 🤔 🚩 ❗ ❓ ⁉️ ⛔
# Fragen an Simon, Calls mit Simon, Allgemeine Fragen an Akademiker, PR-Kontakte

-----------
## Fragen an Simon:


- Was war hiermit gmeint?:
    - Um MCMC zu testen, von Hand gucken, welche Transits passen, die anderen 
      aus den daten werfen
    - dann nur einen Param leicht veraendern
    - und nur den fitten

- Kann/solte man die Walker mit schlechter acceptance aussortieren?

- Welche Parameter kann ich bestimmen, wenn ich nur die TT fitte?
  - Bei einem einzelnen Planeten
    - T0, P, a? 
    - Also m aus P, a? 
    - Was mit e? 
    - Was mit Omega, omega, pomega?
    - r, i: Nein!?!
  - In TOI4504?

   
- Hast du Literatur zu multiple Transits / Fallunterscheidung bei gegenseitiger Ueberdeckung von 3 oder mehr Kreisen?


### Erster Call mit Simon, 24.04.25:

### Call mit Simon, 12.05.25:

- Zeigt der gute Plot wirklich beide Transits?
- Die beiden Transits nebeneinander mit der gefitteten Linie plotten. (Siehe 
  cal_phot.py 227-228)
- for_Uli_mcmc.py debuggen, damit es nicht abstuerzt, wenn man nur einen 
  einzigen Parameter fittet (Zeilen 133-141 und 171-179)
- Aus dem CSV alle Zeilen entfernen, die nicht nahe an den beiden Transits 
  liegen.
- Den Flux der beiden Transits getrennt voneinander auf 1 normieren.
- Zunaechst nur einen Parameter fitten.
  - Reihenfolge:
      - P
      - T0
      - Radiusverhaeltnis R-Planet / R-Stern
      - inc
      - limb darkening

### Next Steps nach Call mit Simon, 04.06.25:

- Mail an Simon schreiben:
    - Github pfad zum Clonen
    - Relevante Verzeichnisse erklaeren

- Simon zum aendern auf GitHub berechtigen?

- Uncertainties korrekt skalieren
- mit mcmc spielen
- Simon schickt Code Update fuer priors

### Next Steps nach Call mit Simon, 16.06.25:

- Zuerst nach d-Transits in den alten Daten suchen:
    - BLS benutzen! https://github.
  com/spacetelescope/tessworkshop_tutorials/blob/master/bls/bls-tutorial.ipynb
    - (Ganz eventuell auch Lightkurve benutzen: https://lightkurve.github.
      io/lightkurve/tutorials/1-getting-started/what-are-periodogram-objects.html )
    - b- und c-Transits aus den Daten entweder ganz entfernen oder rausrechnen.
    - Letzteres haette den Vorteil, dass dann nicht potenziell d-Transits als 
      Kollateralschaden mit entfernt werden.
    - Moegliche Tests:
        - Finde ich b-Transits in den gesamten Daten?
        - Finde ich die beiden d-Transits in den neuen Daten?
        - Finde ich Transits in von mir simulierten Daten?

- CurveSimulator mit MCMC verknuepfen:
    - In Funktion log_likelihood die Berechnung von 
  residuals_phot_sum_squared durch meinen Code ersetzen
    - In "sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, 
      args=(phot_data, spec_data, para, fitting_indices, transformer)"
        enthaelt args die parameter von log_probability (ausser theta = Vektor der Parameter, die MCMC veraendern darf)
        reichen als args data (TESS Messwerte) und para (Input, also Keplerparameter etc)
    - Meine simulierten Fluxdaten muss ich so anpassen, dass sie in jedem 
      Sektor ueber den gleichen Zeitraum integrert sind, wie die jeweiligen TESS-Daten
    - Simon schickt Update von seinem Code


### Call mit Simon, 03.07.25:

- Am besten alle Flux-Daten ausserhalb von Transits vom MCMC fernhalten.
- Waere gut, wenn Curvesimulator zwischen den Transits moeglichst schnell 
  integrieren wuerde. (Veraenderliches dt, oder zumindest Vielfache von dt).
- Residuen durch Uncertainty teilen, damit die Residuen von ungenauen Daten 
  weniger gewichtet werden.
- Next Date: spontan/kurzfristig, falls ich feststecke.


### Call mit Simon, 21.08.25:

- Alle Messungen auf 1800s binnen

- Ulis Ideen fuer gut befunden:
    - Zuerst Check: veraendern sich c und d TT, wenn ich b rausnehme?
    - Planet b: Alle TESS-Daten (ausser den c- und d-Transits) nehmen und damit b fitten.
    - probability zusaetzlich (und am Anfang ausschliesslich) aus Abstand zu den gemessenen TT berechnen
    - Planet e: gibt es bestimmt, aber puh.
    - RV, allerdings sind die Messdaten schlecht  

     
- Erstmal statt MCMC mit einfacherer Methode nur die TT fitten:
    - Levenberg Marquardt LMfit
        - https://pypi.org/project/lmfit/
        - https://lmfit.github.io/lmfit-py/

- Beispiele fuer Systeme mit TTVs, die ich nutzen kann, um CureSimulator zu 
testen:
    - Am besten KEPLER-9. Viele Paper. Klassiker.
    - TOI-1130: https://ui.adsabs.harvard.edu/abs/2023A%26A...675A.115K/abstract
        - Allerdings stammen nur ein Teil der benutzten Photometriedaten von 
      TESS :(
        - In den Referenzen dieses Papers nach geeigneteren Systemen suchen?
    - TOI-1136

- Wenn ich nur Photometrie habe und keine Planeten mit Wechselwirkungen, kann 
ich nur diese Parameter bestimmen:
    - T0, P, a, i, r (weil Stern-r fix)
    - Das ist z.B. relevant, wenn ich nur Planet b alleine fitte
    - P und a separat, weil die Masse von b unbekannt
    - eigentlich a / r-Stern und r / r-Stern statt a und r, aber weil ich ja 
      alle Sternparameter fix halte, passt das so.

- MCMC testen:
    - von hand gucken, welche Transits passen, die anderen aus den daten 
  werfen, dann nur einen Param leicht veraendern und nur den fitten


### Call mit Simon, 01.10.25:

#### Was ich gemacht habe
- gezeigt, dass man b nicht weglassen darf (oder zumindest nicht die Masse 
  von b)
- b erfolgreich gefittet
- Sektor 95: 1c, 1d Transit, habe jetzt 13+4.
- CurveSimulator kann jetzt mit MCMC auch TT statt Flux fitten.
- CurveSimulator kann jetzt auch  mit LMfit die TT fitten.
- Ich habe KEPLER9 problemlos gefittet
- RV hab ich noch ignoriert

#### Simons Empfehlungen fuer die naechsten Schritte
- b weglassen, dafuer zusaetzlich zu c und d auch die Sonnenmasse mit prior 
fitten
- Vitkova plot (mit den roten Punkten) der Verspaetungen der TT gegenueber 
  einer konstanten Periode reproduzieren
- TT mit Vitkova-Parametern reproduzieren
  - Nur die ersten 11 TT von c nehmen
  - Vitkova Parameter nehmen
  - Eine einzige Simulation machen
  - Deltas angucken
- Mit MCMC Vitkova Parameter reproduzieren
  - Nur die ersten 11 TT von c nehmen
  - Vitkova Parameter als Startwerte nehmen
  - MCMC mit sehr kleiner Streuung um diese Startwerte starten
- Als Teil der MCMC-Results die Deltas zwischen observed und computed TT 
  mit den MaxL Parametern plotten


### Call mit Simon, 13.11.25:
#### Was ich gemacht habe
- Passende Parameter fuer TOI4504, TOI4504c, TOI4504d gefunden
- Konnte Vitkova Prognose nicht reproduzieren
- Konnte Vitkova TTV (Figure 6) nicht reproduzieren

#### Fragen/Antworten
- Kannst du Vitkova-Parameter oder Uli-Parameter mit irgendeiner anderen
  (n-body-)Software verifizieren?
  - Nein

#### Next Steps:
- MCMC mit Streuung starten, die der erwarteten Streuung entspricht
- Weitere andere System fitten, um meine Software zu validieren
- Meine Software testen mit einem dieser Systeme/Paper 
  - Kepler18
    - https://ui.adsabs.harvard.edu/abs/2011ApJS..197....7C/abstract
  - TOI1130
    - https://ui.adsabs.harvard.edu/abs/2023A%26A...675A.115K/abstract
  - TOI-216
    - https://iopscience.iop.org/article/10.3847/1538-3881/ab24ba/pdf
  - Diverse Kepler
    - https://arxiv.org/abs/1201.5415
- Chi Square ausrechnen und bei Fits in die Results schreiben
  - Dann Plots machen, wo ich alle bis auf einen Parameter auf den optimalen 
    Werten festhalte, den freien Parameter auf der x-Achse und Chi Squared 
    auf der y-Achse plotten
  - Siehe Blatt Chi Squared in TOI4504.xlsx
- RV in TOI 4504 (zusaetzlich) mit fitten
- Sternmasse als Prior einbauen
  - Also als zusaetzliche Beobachtung inklusive Standardabweichung
  - Implementierung siehe Code-Schnipsel in Mail von Simon vom 13.11.25
- Simon: verstehen, was mean osculating period ist
- Simon erste Dezemberhaelfte im Ausland und 14.1.-23.1. Pruefungen

### Call mit Simon, 22.12.25:
#### Was ich gemacht habe
- Transits Sektor 97 ✅ 
- MCMC mit Streuung starten, die der erwarteten Streuung entspricht ✅ 
- Chi Square ausrechnen und bei Fits in die Results schreiben ✅
  - Dann Plots machen, wo ich alle bis auf einen Parameter auf den optimalen 
    Werten festhalte, den freien Parameter auf der x-Achse und Chi Squared 
    auf der y-Achse plotten 🤔
- RV in TOI 4504 fitten ✅❗
- RV in TOI 4504 zusaetzlich mit fitten 🤔
- Weitere andere System fitten, um meine Software zu validieren 🤔
- Sternmasse als Prior einbauen 🤔
- Was ist mean osculating period? ✅


#### Fragen/Antworten
- KEPLER9 Paper
  - Section 3: "Finally, for each individual planet we correct the output 
    time by the light-travel time effect."
  - Sollte ich das auch machen? Das heisst, wenn der Planet p 1 Lichtminute 
    Abstand vom Stern hat und Planet q 5 Lichtminuten Abstand vom Stern hat: 
    - _subtrahiere_ ich 1 Minute vom p-Transitszeitpunkt 
    - und subtrahiere 5 Minuten vom q-Transitzeitpunkt?❓

- Kommt aus MCMC auch die Standardabweichung der gefitteten Parameter raus?
    - Ich hab einfach mal mean und std von den samples genommen.
    - War das korrekt oder wie kriegt man die raus? ❓
    - Varianz mit der Autokorrelation multiplizieren?
    - Varianz durch die Anzahl walker dividieren?
- 
- Wie funktionieren Jacobi Koordinaten? ❓

- Was sind die naechsten Schritte? ❓


#### Next Steps:



## Allgemeine Fragen an Akademiker:

- How do I get my project used by others?
    - Publish as a paper?
    - How do I get my project published as a paper?
- Where do astronomers gather? Is https://astronomy.stackexchange.com/ the best place to ask questions?
- Is there a community of astronomers specialized in exoplanets that would be welcoming to a dabbling layman like me?
- How do I reliably find out about new papers with interesting star systems? Searching on "arxiv.org"? Or is there a better way?



## PR-Kontakte:

### Bekannte:
- https://www.noraeisner.com/
- Aiden
- Brian
- Peters Kumpel

### Universitäten:
- CZ:
    - michaela.vitkova@asu.cas.cz, vitkova@mail.asu.cas.cz
    - https://stel.asu.cas.cz/en/extrasolar-planet-research/extrasolar-members/
- Niels Bohr:
    - Uffe
    - https://nbi.ku.dk/english/research/astrophysics/exoplanets/
    - http://www.mindstep-science.org/about_us.html ???
    - https://www.iau.org/administration/membership/individual/11759/ ???
- Aarhus:
    - https://phys.au.dk/sac/our-research/extrasolar-planets/
    - Simon Albrecht, Associate Professor, Department of Physics and Astronomy, +4587155702, albrecht@phys.au.dk
        - https://pure.au.dk/portal/en/persons/albrecht%40phys.au.dk
        - https://space.mit.edu/home/albrecht/
    - Professor Emeritus Jørgen Christensen-Dalsgaard
        - Stellar Astrophysics Centre (SAC) an der Aarhus Universität
- SDU (Odense):
    - https://portal.findresearcher.sdu.dk/da/persons/toby-hinse, +4565507131, toch@cp3.sdu.dk
- DTU:
    - http://www.exoplanets.dk/
    - https://orbit.dtu.dk/en/persons/lars-christian-astrup-buchhave, buchhave@space.dtu.dk
    - Office: +45 4525 9500, office@space.dtu.dk
    - https://www.space.dtu.dk/afdelinger/astrofysik-og-atmosfaerens-fysik


### Organisationen:
- https://www.aavso.org/
- https://www.astropy.org/

### Foren
- https://astronomy.stackexchange.com/
- https://space.stackexchange.com/
- https://forum.astronomisk.dk/
- https://www.astrobin.com/forum/
- https://forum.astronomie.de/

- https://www.reddit.com/r/exoplanets/
- https://www.reddit.com/r/TheExoplanetsChannel/
- https://www.reddit.com/r/askastronomy/
- https://www.reddit.com/r/Astronomy/
- https://www.reddit.com/r/astrophysics/
- https://www.reddit.com/r/space/
