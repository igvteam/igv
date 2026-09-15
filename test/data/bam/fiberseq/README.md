# Fiber-seq test BAMs

Copied from the fibertools-rs repository (https://github.com/fiberseq/fibertools-rs, `tests/data/`), MIT licensed.
Indexes were generated with `samtools index`. All files are aligned to hg38.

| File | Tags | Contents |
|---|---|---|
| `nuc_example.bam` | legacy `ns`/`nl`, `as`/`al`/`aq` | 1 read, chr1:3,367,340-3,371,137, MSPs with FIRE qualities |
| `msp_nuc.bam` | legacy `ns`/`nl`, `as`/`al` | 1 read, chr4:3,074,877-3,082,324 |
| `ma_spelled.bam` | `Ma` | The same read as `msp_nuc.bam`, encoded with molecular annotation tags |
| `NAPA.bam` | `MA`/`AQ` (fibertools 0.10-0.12 spelling), `fire` qualities | 154 reads, chr19:47,480,180-47,545,185 |
