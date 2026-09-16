Synthetic RNA-seq BAM for Sashimi intron compression testing
============================================================

egfr.bam / egfr.bam.bai   Synthetic reads over EGFR (hg19), NM_005228 exon structure from UCSC refGene.
                          27 introns: longest 122,920 bp (intron 1), median 1,452 bp, shortest 99 bp.
                          Exon body coverage (100bp reads, depth ~2 per position) plus one junction per
                          intron with depth 20..46 so the arcs carry distinct labels.
                          All reads are plus strand with XS:A:+ and sequence of all A's -- coverage is flat
                          by construction, so this is for testing layout, not for looking at real data.
                          Genome: hg19.   Locus: chr7:55,086,000-55,280,000  (~193 kb)

                          NOTE: the gene is longer than the default alignment visibility window (30 kb).
                          Raise "Visibility range threshold (kb)" in Preferences > Alignments (or the RNA
                          tab) to ~300 before the reads and junctions will load.

                          To regenerate: take the NM_005228 exonStarts/exonEnds from UCSC refGene, emit
                          100 bp reads tiled every 10 bp across each exon (depth 2), plus one read per
                          intron with CIGAR 50M<intron length>N50M at depth 20 + intron index, then sort
                          and index with samtools.

heart.bam / .bai          The real RNA-seq test BAM (hg19, SLC25A3, chr12:98,986,000-98,998,500) from
                          igvteam/igv-data.  All its introns are under 2 kb, so intron compression leaves
                          it essentially unchanged -- useful as the "short introns are untouched" case.
