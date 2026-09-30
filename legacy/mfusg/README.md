# Raven - MODFLOW-USG coupling (L. Scantlebury), preserved unchanged

These twelve files are the MODFLOW-USG groundwater coupling as distributed with Raven 4.15, byte for byte. They are kept
here, outside the build, because the MODFLOW 6 coupling now provides `CGroundwaterModel` and the `.rvg` parser in the main
source folder, and both implementations cannot be compiled into one program under the same class name.

The USG coupling was never active in a standard build: it needs the MODFLOW-USG library and the `_MODFLOW_USG_` switch in
`RavenInclude.h`. To use it, build Raven 4.15 (or earlier) as released.

Its central idea - linking the top layer of an existing model to Raven HRUs through a table of overlap weights
(`:OverlapWeights`, [HRU ID] [cell] [w = A_overlap / A_cell]) - lives on in the MODFLOW 6 coupling: the existing-model mode
reads the same table, and answers the open question noted in the parser ("should cell weights be forced to add to 1? Or
HRU?") by normalising the overlaps of each HRU, so that each HRU's recharge volume reaches MODFLOW exactly.
