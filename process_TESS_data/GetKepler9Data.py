# Code von Perplexity zum Downnload von Kepler-9 Fluxdaten von Kepler

import lightkurve as lk
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def download_kepler9_flux_data():
    """
    Lädt alle verfügbaren Kepler Flux-Daten für Kepler-9 herunter

    Kepler-9 Identifikatoren:
    - KIC: 3323887
    - KOI: 377
    - Kepler ID: kplr003323887
    - Koordinaten: RA 19h 02m 17.76s, Dec +38° 24′ 03.2″
    """

    print("=== KEPLER-9 FLUX DATEN DOWNLOAD ===")
    print("Target: Kepler-9 (KIC 3323887)")

    # 1. Suche nach allen verfügbaren Kepler Beobachtungen
    print("\n1. Suche nach Kepler Beobachtungen...")
    search_result = lk.search_lightcurve("Kepler-9", mission="Kepler")
    print(f"   Gefunden: {len(search_result)} Kepler Beobachtungen")

    # Zeige Details der Beobachtungen
    print("\n2. Beobachtungsdetails:")
    for i, obs in enumerate(search_result):
        print(f"   {i+1:2d}. {obs.mission} - {obs.year} - Cadence: {obs.exptime}s")

    # 3. Download aller Lichtkurven
    print("\n3. Lade alle Lichtkurven herunter...")
    print("   HINWEIS: Dies kann mehrere Minuten dauern!")

    lc_collection = search_result.download_all()
    print(f"   Download abgeschlossen: {len(lc_collection)} Lichtkurven")

    # 4. Verarbeitung und Speicherung
    print("\n4. Verarbeite und speichere Daten...")

    for i, lc in enumerate(lc_collection):
        if lc is not None and len(lc) > 0:
            # Bereinige NaN-Werte
            lc_clean = lc.remove_nans()

            # Konvertiere zu DataFrame
            df = lc_clean.to_pandas()
            df['quarter'] = i + 1

            # Speichere einzelnes Quarter
            filename = f"kepler9_quarter_{i+1:02d}.csv"
            df.to_csv(filename, index=False)

            print(f"   ✓ {filename}")
            print(f"     Zeit: {lc_clean.time.min():.2f} - {lc_clean.time.max():.2f} BKJD")
            print(f"     Punkte: {len(lc_clean):,}")
            print(f"     Flux: {lc_clean.flux.min():.0f} - {lc_clean.flux.max():.0f} e-/s")

    # 5. Kombiniere alle Quarters
    print("\n5. Kombiniere alle Quarters...")
    combined_lc = lc_collection.stitch()
    combined_lc_clean = combined_lc.remove_nans()

    # Speichere kombinierte Daten
    combined_df = combined_lc_clean.to_pandas()
    combined_df.to_csv("kepler9_all_data_combined.csv", index=False)

    print("   ✓ kepler9_all_data_combined.csv")
    print(f"     Gesamtzeitspanne: {combined_lc_clean.time.min():.2f} - {combined_lc_clean.time.max():.2f} BKJD")
    print(f"     Gesamtdauer: {(combined_lc_clean.time.max() - combined_lc_clean.time.min()).value:.1f} Tage")
    print(f"     Gesamtpunkte: {len(combined_lc_clean):,}")

    return combined_lc_clean, lc_collection

# Ausführung
if __name__ == "__main__":
    kepler9_lc, quarters = download_kepler9_flux_data()
