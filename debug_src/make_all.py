def make_all(inputs, output):
    """
    Args:
        inputs: list of filenames
        output: string (name of output file)

    The function reads all input files and concatenates them (inputs[0] first).
    All lines that import code (starting with "import" or "from")
    are moved to the top and duplicates are removed.
    The result is stored in output.
    """
    import_lines = []
    code_lines = []
    seen_imports = set()

    for fname in inputs:
        with open(fname, "r", encoding="utf-8") as f:
            for line in f:
                stripped = line.lstrip()
                if stripped.startswith("import ") or stripped.startswith("from "):
                    # Use the stripped content as key to detect duplicates
                    if stripped not in seen_imports:
                        seen_imports.add(stripped)
                        import_lines.append(stripped.rstrip("\n"))
                else:
                    code_lines.append(line.rstrip("\n"))

    import_lines = [line for line in import_lines if not line.startswith("from .cs_")]
    import_lines = [line for line in import_lines if not line.startswith("from curvesimulator")]
    import_lines = [line for line in import_lines if not line.startswith("import fcntl")]

    with open(output, "w", encoding="utf-8") as out:
        for line in import_lines:
            out.write(line + "\n")
        out.write("\n")
        for line in code_lines:
            out.write(line + "\n")



path = "../src/curvesimulator/"
inputs = [path + "cs_animation.py",
          path + "cs_bodies.py",
          path + "cs_body.py",
          # path + "cs_flux_data.py",
          path + "cs_lightcurve.py",
          path + "cs_manual_fit.py",
          path + "cs_mcmc.py",
          path + "cs_parameters.py",
          path + "cs_physics.py",
          path + "cs_rebound.py",
          path + "cs_results.py",
          path + "curvesim.py",
          "../src/run_curvesim.py",
          ]

make_all(inputs, "../src/all.py")