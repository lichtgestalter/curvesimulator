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
    - starten alle mit max_delta=0.0047, mean_delta=0.0017
    - verbessern sich danach nicht oder nur unwesentlich
    - nelder noch am besten mit max_delta=0.0028, mean_delta=0.
      0015
    - die Deltas liegen zwar deutlich unter den von mir angegebenen 
      Uncertainties der TT, 
    - aber warum gehen die Deltas nicht noch naeher an Null?
    - Ich habe mal Uncertainties der TT um Faktor 1000 verkleinert. Das 
      hat nichts geaendert.
  - Method differential_evolution:
    - max_delta=0.0047   mean_delta=0.0017 nach 1. Iteration
    - max_delta=0.1190   mean_delta=0.0565 nach 5k Iterationen
    - max_delta=0.0452   mean_delta=0.0193 nach 20k Iterationen
 
        
- Konvergiert LMFit wenn Startwerte = Richtige Werte + cirka 1%?
  - `TOI-4504_simX024_01.ini`
  - Methode nelder findet den ersten TT (direkt einen Tag nach StartDate) nicht
  - Methode powell findet den ersten TT nach ca. 60 Iterationen, aber kommt 
    trotzdem nicht nahe an die Parameter, mit denen die Daten erzeugt 
    wurden.
      - max_delta=1.4991, mean_delta=0.4837 [days] 

- Konvergiert LMFit wenn ich mehr Parameter festhalte?
  - `TOI-4504_simX024_03.ini`
  - d.mass fitten, alle anderen Parameter auf 
    richtigen Wert festgesetzt
    - konvergiert _fuer jeden Startwert_ sehr schnell genauso nah an die 
      richtigen Werte, als 
      haette ich wie oben 9 Parameter frei und alle freien Parameter starten 
      mit ihrem richtigen Wert.


- Konvergiert LMFit fuer ein 1-Planeten-System?


- Konvergiert LMFit fuer b (fest) und c (fitten), aber ohne d?



### MCMC fit von Flux inklusive Sektor 95
...

### MCMC fit von KEPLER-9
...

### Simulierten Flux mit LMfit fitten
- erfordert einige Programmier-Arbeit
 
### ab vergessen, as ich hier wollte
- X024 Daten fitten
- X024 ist ein MCMC fit mit 11 Parametern auf die TT bis Sektor 94 (12c, 3d)
  - d: m e P O o ma
  - c: m e P   o ma

