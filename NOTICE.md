# Copyright and licensing

## Raven–MODFLOW 6 groundwater coupling
Copyright (c) 2026 Rezgar Arabzadeh. All rights not granted by the licenses below are reserved.

| Part | Copyright | License |
|---|---|---|
| Raven source code (unchanged parts) | the Raven Development Team | Artistic License 2.0 (`LICENSE`) |
| Coupling source code: `GroundwaterModel.*`, `GWExternal.cpp`, `GWGeometry.*`, `GWGrid.*`, `GWUtil.h`, `MF6Engine.*`, `ParseGWFile.cpp` | Rezgar Arabzadeh | Artistic License 2.0 (`LICENSE`) |
| Changes to existing Raven files (marked "Modified 2026 by Rezgar Arabzadeh") | the Raven Development Team; changes by Rezgar Arabzadeh | Artistic License 2.0 (`LICENSE`) |
| Test and audit scripts (`tests/`, `examples/Liard_groundwater/tools_*.py`) | Rezgar Arabzadeh | Artistic License 2.0 (`LICENSE`) |
| Manual and user guide (`docs/`) | Rezgar Arabzadeh | Creative Commons Attribution 4.0 International (`docs/LICENSE-docs.md`) |
| Liard River model inputs (`examples/Liard_groundwater/`) | their original authors | not relicensed; see `examples/Liard_groundwater/DATA_NOTICE.md` |
| MODFLOW-USG coupling (`legacy/mfusg/`, unchanged from Raven 4.15) | the Raven Development Team (L. Scantlebury) | Artistic License 2.0 (`LICENSE`) |
| Synthetic MODFLOW 6 example model (`examples/Liard_existing_model/model/`) | Rezgar Arabzadeh | Artistic License 2.0 (`LICENSE`) |
| MODFLOW 6 | U.S. Geological Survey | not included; loaded at run time from the user's own installation |

## Differences from the Standard Version of Raven
As the Artistic License 2.0 asks for modified versions, the differences from Raven 4.15 are documented:
file by file in the manual (chapter "Developer guide and limitations", "Changes to existing Raven files") and in full in
the commit history of this branch (`git diff` against Raven's `main`). The modifications are
offered to the Raven Development Team for inclusion in Raven under the Artistic License 2.0.

## Third-party names
Raven and MODFLOW are the names of their respective projects; their use here describes compatibility only.
