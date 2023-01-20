package org.broad.igv.sam;

import htsjdk.samtools.Cigar;
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
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class SupplementaryAlignment implements Locatable {
    private static final Logger log = LogManager.getLogger(SupplementaryAlignment.class);

    public final String chr;
    public final int start;

    public Strand getStrand() {
        return strand;
    }

    private final Strand strand;
    private final Cigar cigar;

    public final int mapQ;
    public final int numMismatches;

    //potentially expensive with very long reads, compute it lazily
    private Integer lenOnRef = null;

    public SupplementaryAlignment(String chr, int start, Strand strand, Cigar cigar, int mapQ, int numMismatches){
        this.chr = chr;
        this.start = start;
        this.strand = strand;
        this.cigar = cigar;
        this.mapQ = mapQ;
        this.numMismatches = numMismatches;
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
        return chr + ":" + Globals.DECIMAL_FORMAT.format(start + 1) + "-" + Globals.DECIMAL_FORMAT.format(start + getLenOnRef())
                + " (" + strand.toShortString() + ") = " + Globals.DECIMAL_FORMAT.format(getLenOnRef()) + "bp  @MAPQ " + mapQ + " NM" + numMismatches;
    }

    public static List<SupplementaryAlignment> parseFromSATag(String saTag){
        return Arrays.stream(Globals.semicolonPattern.split(saTag))
                .map(SupplementaryAlignment::fromSingleSaTagRecord)
                .sorted(LEADING_CLIP_COMPARATOR)
                .collect(Collectors.toList());
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
        switch(strand){
            case NONE:  // if there's no strand it's not aligned, and we can assume it's in cycle order
            case POSITIVE:
                return counts.getLeft();
            case NEGATIVE:
                return counts.getRight();
            default: throw new IllegalStateException(strand + " is not an expected value for strand");
        }
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
        return start + getLenOnRef();
    }

    public int getLenOnRef() {
        if(lenOnRef == null) {
            lenOnRef = cigar.getReferenceLength();
        }
        return cigar.getReferenceLength();
    }

    static class SupplementaryNeighbors {
        final Alignment alignment;
        final SupplementaryAlignment previous;
        final SupplementaryAlignment next;
        public SupplementaryNeighbors(Alignment alignment,
                                      SupplementaryAlignment previous, SupplementaryAlignment next) {
            this.alignment = alignment;
            if(alignment.isNegativeStrand()){
                this.next = previous;
                this.previous = next;
            } else {
                this.next = next;
                this.previous = previous;
            }
        }
        public SupplementaryAlignment previousIgnoreOrientation (){
            return alignment.isNegativeStrand() ? next : previous;
        }

        public SupplementaryAlignment nextIgnoreOrientation (){
            return alignment.isNegativeStrand() ? previous : next;
        }

    }

    public static class SupplementaryGroup {
        private final Alignment alignment;
        private final TreeSet<SupplementaryAlignment> readOrder;
        private final TreeSet<SupplementaryAlignment> positionOrder;

        private final Adapter adapter;

        public SupplementaryGroup(Alignment alignment){
            final List<SupplementaryAlignment> supplementaryAlignments;
            if( alignment instanceof SAMAlignment){
                supplementaryAlignments = ((SAMAlignment) alignment).getSupplementaryAlignments();
            } else {
                final Object rawSATag = alignment.getAttribute(SAMTag.SA.name());
                supplementaryAlignments = rawSATag == null ? null : new ArrayList<>(parseFromSATag(rawSATag.toString()));
            }

            this.alignment = alignment;
            this.adapter = new Adapter(alignment);
            final List<SupplementaryAlignment> combined = new ArrayList<>(supplementaryAlignments);
            combined.add(adapter);
            readOrder = new TreeSet<>(LEADING_CLIP_COMPARATOR);
            readOrder.addAll(combined);
            positionOrder = new TreeSet<>(SortOption.POSITION_COMPARATOR);
            positionOrder.addAll(combined);
        }

        public SupplementaryAlignment getNextInRead(SupplementaryAlignment alignment){
            return readOrder.higher(alignment);
        }

        public SupplementaryAlignment getPreviousInRead(SupplementaryAlignment alignment){
            return readOrder.lower(alignment);
        }
        public SupplementaryAlignment getNextPosition(SupplementaryAlignment alignment){
            return positionOrder.higher(alignment);
        }

        public SupplementaryAlignment getPreviousPosition(SupplementaryAlignment alignment){
            return positionOrder.lower(alignment);
        }

        public Adapter getAdapter() {
            return adapter;
        }

        public Iterator<SupplementaryAlignment> iterateInReadOrder(){
            return readOrder.iterator();
        }

        public Iterator<SupplementaryAlignment> iterateInPositionOrder(){
            return positionOrder.iterator();
        }

        public Stream<SupplementaryAlignment> streamInReadOrder(){
            return readOrder.stream();
        }
        public Stream<SupplementaryAlignment> streamInPositionOrder(){
            return positionOrder.stream();
        }

        public int size(){
            return readOrder.size();
        }

        private static class Adapter extends SupplementaryAlignment{
            public Adapter(Alignment a) {
                super(a.getChr(), a.getStart(), a.getReadStrand(), a.getCigar(), a.getMappingQuality(), a.getAttribute(SAMTag.NM.name()) == null ? 0 : (int)a.getAttribute(SAMTag.NM.name()));
            }
        }
    }
}