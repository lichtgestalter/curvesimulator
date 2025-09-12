- Simulierte Daten fitten lassen
  - Generiere mit den X024ML Parametern Simulierten Flux und TT Sim002
      - TOI4504c Inklination von 89.69 auf 89.79, damit Planet c nicht zu 
        frueh rausrutscht (also aufhoert Transits zu erzeugen)
      - Configfile TOI-4504_X024_maxl.ini 
      - TOI-4504_X024_maxl_dt20000.json erzeugt
      - Mit results2list() und df2csv() in result_evaluation.py erzeugt:
          - TT_TOI4504d.csv
          - TT_TOI4504c.csv
      - Daraus die TT rausgesucht, die ungefaehr den gemessenen TT (13c, 4d) 
        entsprechen und diese in ../data/simulations/tt_simX024.csv gespeichert
      - Configfile TOI-4504_simX024_01.ini aus TOI-4504_X024_maxl.ini erzeugt
          - zunaechst die Body-Parameter nicht veraendert. D.h., LMfit 
            bekommt die optimalen Parameter-Werte bereits als Startwerte 
            serviert.
          - Nur 9 Parameter sind zu fitten:
              - d: m e P   ma
              - c: m e P O ma
          - Method powell:
              - max_delta=0.0047   mean_delta=0.0017 nach 1. Iteration
              - max_delta=0.0039   mean_delta=0.0016 nach ca. 100 Iterationen
          - Method nelder:
              - max_delta=0.0047   mean_delta=0.0017 nach 1. Iteration
              - max_delta=0.0028   mean_delta=0.0015 nach ca. 100 Iterationen
          - Method lbfgsb:
              - max_delta=0.0047   mean_delta=0.0017 nach 1. Iteration
              - max_delta=0.0047   mean_delta=0.0017 am Ende
          - Method cg:
              - max_delta=0.0047   mean_delta=0.0017 nach 1. Iteration
              - max_delta=0.0047   mean_delta=0.0017 nach ca. 200 Iterationen
          - Method cobyla:
              - max_delta=0.0047   mean_delta=0.0017 nach 1. Iteration
              - max_delta=0.0049   mean_delta=0.0017 am Ende
          - Method bfgs:
              - max_delta=0.0047   mean_delta=0.0017 nach 1. Iteration
              - max_delta=0.0047   mean_delta=0.0017 nach ca. 200 Iterationen
          - Method tnc:
              - max_delta=0.0047   mean_delta=0.0017 nach 1. Iteration
              - max_delta=0.0047   mean_delta=0.0017 nach ca. 200 Iterationen
          - Method trust-constr:
              - max_delta=0.0047   mean_delta=0.0017 nach 1. Iteration
              - max_delta=0.0047   mean_delta=0.0017 am Ende
          - Method differential_evolution:
              - max_delta=0.0047   mean_delta=0.0017 nach 1. Iteration
              - max_delta=0.1190   mean_delta=0.0565 nach ca. 5200 Iterationen
 
 
        
- Ablage unter data/simulations
- Konvergiert LMFit wenn Startwerte = Richtige Werte? 
- Konvergiert LMFit wenn Startwerte = Richtige Werte + x% Rauschen? 
- Verschiedene Methoden ausprobieren (Nelder, ...)
- LMfit mit Fluxdaten laufen lassen? (erfordert neuen Code)
 