**Copy number color scales and expected ploidy**

***Background***

Written 2026-11-24 while adding "color genotypes by FORMAT field" to variant tracks (branch
`color-genotypes-by-format`).  That work gave the FORMAT field CN a default scale matching segmented copy number
tracks, which exposed a long standing flaw shared by both: the scale's neutral point assumes two copies
everywhere.  Nothing was changed for this.  This note records the problem so it can be picked up later.

***The problem***

- Segmented copy number tracks are colored by `IGVPreferences.getColorScale(DataType.COPY_NUMBER)`: blue --
  neutral -- red, neutral from -0.1 to 0.1 and saturated at +/- 1.5 (`getDefaultColorScale`).  The values are
  taken to be log2(copy number / 2), so neutral means two copies.
- The CN default in `VariantTrackMenuHelper.copyNumberScale()` uses the same colors on an absolute copy number:
  neutral at 2, full blue at 0, full red at 6 and above.
- Where the expected copy number is not 2, the normal state is drawn as a loss.  On chrY, and on chrX in a
  male sample, one copy is normal: a seg value of -1 and a CN of 1 both draw clearly blue, and a true
  duplication to two copies draws neutral.
- The 1000 Genomes phase 3 chrY calls show it directly -- nearly every genotype there has CN 1.  Example loci
  are `chrY:2,892,962-2,900,130` and `chrY:3,248,822-3,253,727` in `test/sessions/json/vcf-session.json` (raise
  the visibility window to see the larger events).

***What a fix would need***

Not designed.  Points to settle if this is taken up:

- The expected ploidy per chromosome, or per region.  For human it depends on the sample's sex, and the
  pseudoautosomal regions of X and Y are diploid in both.  Other organisms differ (haploid genomes, ZW systems,
  polyploids), so it cannot be hardcoded to X and Y.
- Where the sample's sex would come from.  IGV does not know it.  A sample information attribute is one
  possible source; a VCF's own header or a per-track setting are others.
- Whether to shift the scale or the values.  For CN the neutral point could move to the expected copy number.
  For seg data the file's values are already ratios against whatever reference the caller used, which may or
  may not have been sex-matched, so IGV cannot tell whether a correction is needed.
- The CN scale is linear in copy number between its anchor points, where a seg track is logarithmic.  CN 3 is
  25% red where the equivalent seg value (0.58) is about 39%.  A ploidy-aware scale could be made logarithmic
  at the same time: log2(CN / expected).
