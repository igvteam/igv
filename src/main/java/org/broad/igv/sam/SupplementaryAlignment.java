package org.broad.igv.sam;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.Locatable;
import org.broad.igv.Globals;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class SupplementaryAlignment implements Locatable {
    private static final Logger log = LogManager.getLogger(SupplementaryAlignment.class);

    public final String chr;
    public final int start;

    private final Strand strand;
    private final Cigar cigar;

    public final int mapQ;
    public final int numMismatches;

    //potentially expensive with very long reads, compute  lazily
    private Integer lenOnRef = null;  //number of reference bases covered
    private Integer numberOfAlignedBases = null; //number of bases in read which are aligned

    public SupplementaryAlignment(String chr, int start, Strand strand, Cigar cigar, int mapQ, int numMismatches){
        this.chr = chr;
        this.start = start;
        this.strand = strand;
        this.cigar = cigar;
        this.mapQ = mapQ;
        this.numMismatches = numMismatches;
    }

    public Strand getStrand() {
        return strand;
    }

    public static int getInsertionIndex(final Alignment alignment, final List<SupplementaryAlignment> supplementaryAlignments) {
        final SupplementaryAlignment alignmentAsSupplementary = new SupplementaryAlignment(null, 0, alignment.getReadStrand(), alignment.getCigar(),
                0, 0);  //We're just using this as a placeholder for sorting since Alignment
        // and SupplementaryAlignment don't share any common class
        int insertionIndex = Collections.binarySearch(supplementaryAlignments, alignmentAsSupplementary, LEADING_CLIP_COMPARATOR);
        if(insertionIndex > 0) {
            log.warn("Saw 2 SupplementaryReads with the same leading clip length.");
        } else {
            insertionIndex = -1* (insertionIndex + 1);
        }
        return insertionIndex;
    }

    public static SupplementaryNeighbors getAdjacentSupplementaryReads(final Alignment alignment) {
        final List<SupplementaryAlignment> supplementaryAlignments;

        final SupplementaryNeighbors supplementaryRenderingInfo;

        if( alignment instanceof SAMAlignment){
            supplementaryAlignments = ((SAMAlignment) alignment).getSupplementaryAlignments();
        } else {
            final Object rawSATag = alignment.getAttribute(SAMTag.SA.name());
            supplementaryAlignments = rawSATag == null ? null : new ArrayList<>(parseFromSATag(rawSATag.toString()));
        }


        if(supplementaryAlignments != null && !supplementaryAlignments.isEmpty()) {

            int insertionIndex = getInsertionIndex(alignment, supplementaryAlignments);
            if( insertionIndex < 0 ){
                log.warn("Encountered a group of supplementary alignments that couldn't be sorted unambiguously." +
                        "\nSee reads with the name: " + alignment.getReadName());
                return null;
            }
            //Is this the first or last in a set of supplementary reads by position along the original read.
            // ex:
            // --- Match / xxxx Clip
            // original unaligned read:
            // -----------------------------
            // split into 3
            // xxxxxx--------xxxxxxxxxxxxxx : not first or last
            // ------xxxxxxxxxxxxxxxxxxxxxx: first
            // xxxxxxxxxxxxxx--------------: last

            final boolean firstInRead = insertionIndex == 0;
            final boolean lastInRead = insertionIndex == supplementaryAlignments.size();
            supplementaryRenderingInfo = new SupplementaryNeighbors(alignment,
                    firstInRead ? null : supplementaryAlignments.get(insertionIndex - 1),
                    lastInRead ? null : supplementaryAlignments.get(insertionIndex));
            return supplementaryRenderingInfo;
        } else {
            return null;
        }
    }

    public String printString() {
        // chr6:43,143,415-43,149,942 (-) @ MAPQ 60 NM 763
        // be sure to adjust start by + 1 because SATag is 1 based but IGV internal is 0 based
        return chr + ":" + Globals.DECIMAL_FORMAT.format(start + 1) + "-" + Globals.DECIMAL_FORMAT.format(start + getLengthOnReference())
                + " (" + strand.toShortString() + ") = " + Globals.DECIMAL_FORMAT.format(getLengthOnReference()) + "bp  @MAPQ " + mapQ + " NM" + numMismatches;
    }

    public static List<SupplementaryAlignment> parseFromSATag(String saTag){
        return Arrays.stream(Globals.semicolonPattern.split(saTag))
                .map(SupplementaryAlignment::fromSingleSaTagRecord)
                .sorted(LEADING_CLIP_COMPARATOR)
                .toList();
    }

    public static SupplementaryAlignment fromSingleSaTagRecord(String saTagRecord){
        //SA:(reference name, pos, strand, CIGAR, mapQ, NM;)+
        String[] tokens = Globals.commaPattern.split(saTagRecord);

        //We have to canonicalize the chr name
        String chr = tokens[0];
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        chr = genome == null ? chr : genome.getCanonicalChrName(chr);

        int start = Integer.parseInt(tokens[1]) - 1; // this is 1 based so subtract 1!
        Strand strand = Strand.fromString(tokens[2]);
        Cigar cigar = TextCigarCodec.decode(tokens[3]);
        int mapQ = Integer.parseInt(tokens[4]);
        int numMismatches = Integer.parseInt(tokens[5]);
        return new SupplementaryAlignment(chr, start, strand, cigar, mapQ, numMismatches);
    }


    public int getCountOfLeadingClipping(){
        ClippingCounts counts = ClippingCounts.fromCigar(this.cigar);
        return getCountOfLeadingClipping(counts, this.strand);
    }

    public static int getCountOfLeadingClipping(final ClippingCounts counts, final Strand strand) {
        return switch (strand) {  // if there's no strand it's not aligned, and we can assume it's in cycle order
            case NONE, POSITIVE -> counts.getLeft();
            case NEGATIVE -> counts.getRight();
        };
    }

    /**
     * Sort by how the alignments were positioned in the original read.  Sort by the amount of clipping on the strand in read order.
     */
    public static final Comparator<SupplementaryAlignment> LEADING_CLIP_COMPARATOR = Comparator.comparingInt(SupplementaryAlignment::getCountOfLeadingClipping);

    @Override
    public String getContig() {
        return chr;
    }

    @Override
    public int getStart() {
        return start;
    }

    @Override
    public int getEnd() {
        return start + getLengthOnReference();
    }

    @Override
    public int getLengthOnReference() {
        if(lenOnRef == null) {
            lenOnRef = cigar.getReferenceLength();
        }
        return cigar.getReferenceLength();
    }

    /**
     * get the count of non-clipped bases which are in the read
     * this differs from {@link #getLengthOnReference()} because it includes insertions bases but not deletions
     */
    public int getNumberOfAlignedBases(){
        if( numberOfAlignedBases == null) {
            int length = 0;
            for(CigarElement element : cigar) {
                switch (element.getOperator()) {
                    case M, I, EQ, X -> length += element.getLength();
                    default -> {}
                }
            }
            numberOfAlignedBases = length;
        }
        return numberOfAlignedBases;
    }


    record SupplementaryNeighbors(Alignment alignment, SupplementaryAlignment previous, SupplementaryAlignment next) {
            SupplementaryNeighbors(Alignment alignment, SupplementaryAlignment previous, SupplementaryAlignment next) {
                this.alignment = alignment;
                if (alignment.isNegativeStrand()) {
                    this.next = previous;
                    this.previous = next;
                } else {
                    this.next = next;
                    this.previous = previous;
                }
            }

        }

}