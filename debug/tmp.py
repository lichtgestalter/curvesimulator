import sys
if 'curvesimulator.curvesim' in sys.modules:
    print('deleting curvesimulator.curvesim')
    del sys.modules['curvesimulator.curvesim']
else:
    print('curvesimulator.curvesim was not in sys.modules')


