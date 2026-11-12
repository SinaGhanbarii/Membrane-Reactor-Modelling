# Zotero package: Pd-membrane WGS reactor corpus

Built 12 Nov 2026 from the 45-paper review (literature-matrix.csv) plus the 26 references the review says to retrieve.

## Files
- `zotero_import.py`: creates the collection tree and all items through the Zotero web API (items can sit in several sub-collections). Standard library only.
- `corpus_items.json`: the data the script reads (71 items: 45 corpus papers and 26 to-retrieve).
- `corpus_items.ris`: the 45 corpus papers only, in one file, for a no-key import (File > Import). RIS cannot hold collections, so you place items yourself.
- `collections_map.csv`: which collections each corpus paper goes in, with the file name and the DOI printed in that file.

## Run the script
1. Get your numeric user ID and a write-enabled key at https://www.zotero.org/settings/keys.
2. Dry run (no network): `python zotero_import.py`
3. Create everything: `python zotero_import.py --user-id <ID> --api-key <KEY> --apply`

The script first compares its field names with Zotero's own item templates and stops if any differ. It stops if the root collection already exists, so run it once.

## Collections
Root `Pd-membrane WGS reactor (PhD)` with 10 sub-collections: 01 Kinetics and thermodynamics, 02 Membrane permeation parameters, 03 Gas inhibition, 04 Concentration polarisation and mass transfer, 05 Reactor models, 06 Thermal behaviour, 07 Reviews and reports, 08 Reforming and other reactions, 09 To retrieve, 10 Excluded.

## Tags on every item
`id:Rxx` (review ID), `type:EXP/SIM/MIXED/REVIEW`, `tier:Core/Supporting/Low`. Extra tags: `metadata-unverified`, `off-topic`, `to-retrieve`, `preprint-not-peer-reviewed`.

## Cautions
- The collection assignment is a judgement from the review text. Move items freely.
- DOIs are included only where the DOI is printed in the supplied file. Keiski et al. (1996) has none printed.
- Items tagged `metadata-unverified` (R23, R29, R40-R43 and all 26 to-retrieve items) were copied from title-less files or from other papers' reference lists. Check them against the publisher record.
- PDFs are not attached. The `file` column in the CSV tells you which PDF belongs to which item.
- The script was tested against a stand-in for the API, not the live Zotero service.
