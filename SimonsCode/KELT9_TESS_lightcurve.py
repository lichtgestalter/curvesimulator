
import numpy as np
import pandas as pd

import lightkurve as lk
import matplotlib.pyplot as plt

search_result = lk.search_lightcurve('KELT-9', mission='TESS')
search_result

lc = search_result[0].download()

lc.plot()
plt.show()


time_bjd = lc.time.value  
flux = lc.flux.value      
flux_unc = lc.flux_err.value  # Uli
# flux_unc = lc.flux_unc.value


lightcurve_data = pd.DataFrame({
    "BJD": time_bjd,
    "Flux": flux,
    "Flux_unc": flux_unc
})


lightcurve_data.to_csv("data/KELT9_TESS_phot.csv", index=False)

# print("downloaded and saved lightcurve.csv'")