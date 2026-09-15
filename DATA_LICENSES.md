# Data Sources and Licensing

## Structural coordinates — NOT redistributed

Every structure is referenced by **PDB ID** and fetched at run time from
`https://files.rcsb.org/download/<id>.cif`. Coordinates are cached in `.structure_cache/`, which is
git-ignored and marked read-only on download. **No coordinate file is redistributed with this
repository.**

PDB data are released by the wwPDB under CC0 1.0 (public domain dedication). Citation of the
originating depositors' primary publications remains expected practice and will be provided for
every structure discussed in any manuscript.

## Annotations

- **SIFTS** residue-level UniProt↔PDB alignments, obtained via the RCSB Data API
  (`rcsb_polymer_entity_align`, `provenance_source = SIFTS`). SIFTS is produced by PDBe/EMBL-EBI
  and released under CC BY 4.0.
- **RCSB Search and Data APIs** for entry metadata and unobserved-residue features.
  Derived annotation files under `annotations/` and `data_manifest/` record accessions, residue
  ranges and provenance, not coordinates.

## Software

| Tool | Version | Licence |
|---|---|---|
| P2Rank | 2.5.1 | MIT |
| fpocket | source tag 4.2.3 | GPL-3.0 (bundles qhull, which carries its own licence) |
| gemmi | 0.7.5 | MPL-2.0 |
| numpy / scipy | 2.5.2 / 1.18.1 | BSD-3-Clause |

fpocket is GPL-3.0. It is invoked as an external executable and is neither vendored nor linked into
this repository's code, so no copyleft obligation attaches to these scripts. Anyone redistributing a
built fpocket binary alongside this work must comply with GPL-3.0 themselves.
