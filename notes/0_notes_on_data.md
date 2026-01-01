### TOI4504_transits_sm0_3til94.csv
    Made using get_all_c_d_transits()
    12 c transits and 3 d transits in sectors 3 - 94
    Sector 3: TGLC 1800s. No flux error provided. We set normalized flux error = 0.002
    Sector 6, 9, 12: QLP 1800s
    Sector 61: QLP 200. We multiplied flux error by 3 (at transit, sun came almost into TESS' field of view).
    Other sectors: SPOC 120s

### TOI4504_transits_sm0_3til95.csv
    Added Sector 95 manually
    Now 13 c transits and 4 d transits


### TOI4504_no_transits_sm0_27til94.csv
    def combine_all_flux():
        # read all downloaded SPOC data from sectors >= 27
        # remove c and d transits
        dfs_cleaned_from_cd_transits = get_all_c_d_transits(spoc_only=True, no_transits=True)
        dfs_from_other_sectors = get_other_flux_for_b_transits()
        all_df = pd.concat(dfs_cleaned_from_cd_transits + dfs_from_other_sectors, ignore_index=True)
        all_df = all_df.dropna(subset=["flux"]).sort_values(by="time", ascending=True)
        all_df = all_df[(all_df["flux"] >= 0.98) & (all_df["flux"] <= 1.02)]
        df2csv(all_df, path + "TOI4504_no_transits_sm0_27til94.csv")

### TOI4504_b_transits_bin1200_27til94.csv
    df = csv2df(path + "TOI4504_no_transits_sm0_27til94.csv")
    exposure_time = 120
    bin_scope = 1200
    df_binned = bin_flux(df, exposure_time, bin_scope)
    df_extracted = extract_regular_transits(df_binned, 2458400, 2.42614, 0.33, 0.45)
    df2csv(df_extracted, path + f"TOI4504_b_transits_bin{bin_scope}_27til94.csv")

### TOI4504_b_transits_27til94.csv
    df_extracted = extract_regular_transits(df, 2458400, 2.42614, 0.33, 0.45)
    df2csv(df_extracted, path + f"TOI4504_b_transits_27til94.csv")
