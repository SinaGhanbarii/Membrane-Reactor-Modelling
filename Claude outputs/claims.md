# Claims ledger: reactor geometry figure

| Item | Value | Basis | Source |
|---|---|---|---|
| Membrane outer diameter | 9.525 mm (3/8 in) | published template value, in register | cases.base.membrane_OD (Hla 2015 template) |
| Membrane length (bed length) | 0.54 m | published template value, in register | cases.base.membrane_length |
| Shell inner diameter | 25.4 mm (1 in) | assumed (placeholder, shell ID not stated by Hla) | cases.base.shell_ID |
| Bed voidage | 0.40 | assumed | cases.base.bed_voidage |
| Catalyst particle diameter (drawn to scale) | 0.7 mm | register template, Hla 2015 range 0.5-0.85 mm | cases.base.d_p |
| Annulus gap | 7.9 mm | calculated: (25.4 - 9.525) / 2 | - |
| Bed cross-section | 4.35 cm2 | calculated: pi/4 (D_s^2 - D_m^2) | - |
| Membrane area | 162 cm2 | calculated: pi D_m L (Hla 2015 states 150 cm2; the difference is flagged in the design document) | - |
| Annulus volume | 235 cm3 | calculated: A_bed x L | - |
| Axial compression | about 5:1 | drawing (0.54 m at the cross-section scale would be about 2550 px; drawn 530 px) | - |
| Membrane wall and shell wall thickness | not to scale | drawing | - |
| Permeate inside the tube | model assumption | not confirmed for the H2Site rig | design document, section 2 |
