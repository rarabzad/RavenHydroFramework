# Example: Raven coupled to an existing MODFLOW 6 model

`model/` is a MODFLOW 6 model of the Liard coupled area built independently of Raven (UTM zone 10N, 1.5 km cells rotated
15 degrees, four layers, a one-day steady-state period and 24 monthly periods, two wells with seasonal pumping, 239 river
cells with their Raven subbasin ID in an auxiliary variable, a head boundary at the outlet, its own recharge and ET).
Raven's files are those of `../Liard_groundwater` (forcing and observations are read from there).

Run (two years, about four minutes):

    ../../src/Raven.exe Liard -o output/        # the library from tools/get_mf6.py is found by itself
    python3 ../Liard_groundwater/tools_audit.py output

The run works on a copy of the model in `output/mf6/`; `model/` is never modified. `GWModelSummary.txt` describes the
coupling (model, units, how cells were linked, packages taken over and kept, HRU coverage). See the manual, chapter
"Coupling an existing MODFLOW 6 model".

The model was generated with FloPy (`tests/existing/make_ext_model.py`); FloPy is not needed to run the example.
