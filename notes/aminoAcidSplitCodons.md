# Amino acid rendering: split codons with no neighbouring exon

November 2026. Findings from a review of IGV desktop prompted by an igv.js bug report.

## Background

When zoomed in far enough, feature tracks draw the amino acids of coding exons. A codon can be split by an exon
edge, with one or two of its bases in the neighbouring exon. To translate it, the code borrows the missing bases
from that neighbour.

The neighbour is not always there. A feature can cover only part of a protein, so its first or last exon begins or
ends mid-codon with nothing beside it. Examples:

- UCSC's UniProt bigBed tracks (domains, regions and so on), where each feature is a stretch of a protein.
- Transcripts with an incomplete CDS, such as a GFF/GTF transcript whose first CDS has a non-zero phase.

igv.js crashed in this situation ([igvteam/igv.js#2127](https://github.com/igvteam/igv.js/issues/2127), fixed in
[igvteam/igv.js#2128](https://github.com/igvteam/igv.js/pull/2128)): three of the four cases (leading or trailing
split codon, plus or minus strand) dereferenced the missing neighbour. The igv.js code was derived from IGV
desktop, so desktop was checked for the same problem.

## Reading frame convention

`readingFrame` on an exon is the number of bases of its first codon that lie in the *preceding* exon, in the
direction of transcription (0, 1 or 2). This is what `FeatureUtils.computeReadingFrames` computes and what UCSC's
`exonFrames` holds. `GFFCodec` converts a GFF phase to it with `(3 - phase) % 3`.

So an exon with `readingFrame` 2 has one base of its first codon at its own 5' end; with `readingFrame` 1 it has
two. On the plus strand the 5' end is the exon's start and the preceding exon is `prevExon`. On the minus strand
the 5' end is the exon's end and the preceding exon is `nextExon`.

## Where the code is

- `Exon.computeAminoAcidSequence(genome, prevExon, nextExon)` builds the nucleotide string for the exon, padding
  it with bases from the neighbours, and translates it with `AminoAcidManager.getAminoAcidSequence`.
- `AminoAcidManager.getAminoAcids` translates from the first base of the string on the plus strand, and from the
  last base on the minus strand. It assumes the string is already in frame; leftover bases at the far end are
  ignored.
- `IGVFeatureRenderer.labelAminoAcids` passes the neighbouring exons (or null) and draws the result.

This is structured differently from igv.js, which translates whole codons in a loop and handles the two junction
codons separately.

## Findings

### No crash for a missing neighbour

All four cases check the neighbour for null, so desktop does not have the igv.js crash. The case in the igv.js
report (the minus-strand ADAM17 "Disintegrin" feature in
`https://hgdownload.soe.ucsc.edu/gbdb/hg19/uniprot/unipDomain.bb` at `chr2:9637312-9637362`, frame 0 with two
leftover bases and no exon to the left) is handled: the leftover is left untranslated.

### Out-of-frame translation (fixed)

When `readingFrame > 0` and the preceding exon was missing, the padding was skipped but the bases of the split
codon stayed at the front of the string. Translation then started at the exon edge instead of after those bases,
so every amino acid in the exon was out of frame.

The fix, in `Exon.computeAminoAcidSequence`, drops the `3 - readingFrame` bases of the split codon when there is no
exon to complete it (from the start of the string on the plus strand, from the end on the minus strand) and moves
`readStart` or `readEnd` to match. The split codon is left untranslated, as in igv.js. A non-coding (UTR) preceding
exon is treated the same as a missing one, since its bases are not part of the codon.

`IGVFeatureRenderer.labelAminoAcids` fills the bases of an untranslated split codon, at either edge of the exon and
on either strand, with the alternating amino acid background and no letter. The codon boundaries stay visible even
though the codon is not known.

Regression tests: `ExonTest.testSplitCodonNoPrevExonPos` and `ExonTest.testSplitCodonNoNextExonNeg`. Each compares
an exon with a split leading codon and no neighbour against an exon with frame 0 that starts just past the split
codon; both must give the same amino acids at the same position. `testSplitCodonNonCodingPrevExonPos` and
`testSplitCodonNonCodingNextExonNeg` do the same with a non-coding neighbour in place of the missing one.

### Not fixed

- **Null neighbour sequence.** Three of the four cases use `d.length` without checking that `genome.getSequence`
  returned something for the neighbour's bases; only the plus-strand next-exon case has `if (d != null)`. A null
  there would throw a `NullPointerException`. It is unlikely in practice, since the exon's own sequence has just
  loaded from the same genome.
- **Non-coding neighbour for a trailing codon, minus strand.** The minus-strand previous-exon case does not check
  `isNonCoding()` and would borrow bases from a UTR exon to complete the last codon. The other three cases do
  check. igv.js requires the neighbour to have a reading frame. Whether real data reaches this was not worked out.

### Data caveat: UniProt exonFrames

In the `unipDomain.bb` feature above, `exonFrames` is `0,0,0` for exons of 41, 104 and 122 bases, which cannot be
right for a continuous coding sequence. With those frames each exon is translated from its own edge, so the amino
acids in the downstream exons are out of frame in both igv.js and desktop. Only this one feature was checked. If
the UniProt files do this generally, the amino acids shown on those tracks are unreliable, whatever the code does.
