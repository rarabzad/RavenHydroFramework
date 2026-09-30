# Calibrating the Raven-MODFLOW 6 coupling with Ostrich

`ostIn.txt` calibrates gravel K, gravel Sy and upland-till K against monitoring-well heads (GW_HEAD in
Raven's Diagnostics.csv) with DDS. `Liard.rvp.tpl` is the .rvp with parameter tokens in :AquiferClasses.
Check the line/column of each response against your own Diagnostics.csv before running.

Not run with the Ostrich executable here. The same search was run with the Python DDS driver
`../calib.py` (identical templates, objective and bounds) on the synthetic subbasin-52 test; see
`../calib_log.txt` and the reference document for the results.
