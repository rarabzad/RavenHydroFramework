# Tests of the existing-model mode

Development tests (paths inside the scripts point to the test machine; adjust them to your folders).

- `make_ext_model.py` - builds the synthetic MODFLOW 6 model of the Liard area (FloPy, pyproj; only for building it).
- `run_standalone.py` - runs a MODFLOW 6 simulation through libmf6 alone.
- `setup_case.py` - a Liard Raven case coupled to a model folder.
- `guards.py` - 14 faulty set-ups; each must stop with a clear message.
- `roundtrip.py` - a Raven-generated model re-imported as an existing model must reproduce the generated run.
- `make_weights.py` - :OverlapWeights computed independently (shapely), for the weights test.
- `restart_ext.py` - restarts inside a stress period and before a change of well rates against a continuous run.
- `fuzz_ext.py SEED N` - randomized set-ups (units, link method, stream mapping, seepage, ET, coverage, restarts), checked
  with the audit identities and MODFLOW's own budget file.
- `crstest.cpp`, `crs_check.py` - the built-in projections against PROJ.
