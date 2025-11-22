def rv1():
    from dace_query import Dace
    # Query radial velocities for TOI-4504
    rv_data = Dace.query_radial_velocities('TOI-4504')
    print(f"Found {len(rv_data)} RV measurements")
    print(rv_data[['rjd', 'rv', 'rv_err', 'ins_name']])

def rv2():
    from dace_query.spectroscopy import Spectroscopy

    # Query FEROS observations of TOI-4504
    observed_data = Spectroscopy.query_database(
        limit=1000,
        filters={
            "public": {"is": True},
            "obj_id_catname": {"contains": ["TOI-4504"]},
            "ins_name": {"contains": ["FEROS"]}
        },
        output_format="pandas"
    )

    # Get radial velocity time series
    spectro_ts = Spectroscopy.get_timeseries("TOI-4504", output_format="pandas")

    # Display results
    print(f"FEROS observations: {len(observed_data)}")
    print(spectro_ts[['rjd', 'rv', 'rv_err', 'ins_name']])

    # Optional: Download raw FITS spectra
    Spectroscopy.download(
        "all",
        filters={
            "obj_id_catname": {"contains": ["TOI-4504"]},
            "ins_name": {"contains": ["FEROS"]}
        },
        output_directory="./TOI-4504-data",
        output_filename="TOI-4504-FEROS.tar.gz"
    )

def rv3():
    from astroquery.utils.tap.core import Tap
    # Connect to NASA Exoplanet Archive TAP service
    tap = Tap(url="https://exoplanetarchive.ipac.caltech.edu/TAP")
    # Query for TOI-4504 planets and parameters
    query = """
    SELECT pl_name, pl_orbper, pl_rade, pl_masse
    FROM ps 
    WHERE pl_name LIKE 'TOI-4504%'
    ORDER BY pl_name
    """
    result = tap.search(query)
    results_df = result.to_pandas()
    print(results_df)

def rv4():
    from astroquery.eso import Eso
    # import getpass

    # Initialize ESO query (authentication optional for public data)
    eso = Eso()

    # Query FEROS observations of TOI-4504
    table = eso.query_instrument(
        'FEROS',
        object='TOI-4504',
        start_date='2023-01-01'
    )

    print(f"Found {len(table)} FEROS observations")

    # Download public data (within 1 year of observation, requires authentication for proprietary)
    # data_files = eso.retrieve_data(table['DP.ID'][:10])


# All AI-Code. None of these work.
rv4()