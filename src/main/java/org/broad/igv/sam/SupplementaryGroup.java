package org.broad.igv.sam;

import htsjdk.samtools.SAMTag;
import org.broad.igv.feature.genome.ChromosomeNameComparator;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class SupplementaryGroup {
    private final List<SupplementaryAlignment> readOrder;
    private final List<SupplementaryAlignment> positionOrder;

    private final Adapter primaryRead;
    private final Alignment original;

    public SupplementaryGroup(Alignment alignment) {
        final List<SupplementaryAlignment> supplementaryAlignments;
        if (alignment instanceof SAMAlignment) {
            supplementaryAlignments = ((SAMAlignment) alignment).getSupplementaryAlignments();
        } else {
            final Object rawSATag = alignment.getAttribute(SAMTag.SA.name());
            supplementaryAlignments = rawSATag == null ? null : new ArrayList<>(SupplementaryAlignment.parseFromSATag(rawSATag.toString()));
        }

        this.primaryRead = new Adapter(alignment);
        this.original = alignment;
        final List<SupplementaryAlignment> combined = (supplementaryAlignments == null || supplementaryAlignments.isEmpty())
                ? new ArrayList<>()
                : new ArrayList<>(supplementaryAlignments);

        combined.add(primaryRead);
        readOrder = new ArrayList<>(combined);
        readOrder.sort(SupplementaryAlignment.LEADING_CLIP_COMPARATOR);
        ;

        positionOrder = new ArrayList<>(combined);
        positionOrder.sort(SortOption.POSITION_COMPARATOR);
    }

    /**
     * @return the number of non-clipped bases in the reads.  This includes insertion bases but not deletions.
     */
    public int getBaseCount() {
        return streamInPositionOrder()
                .mapToInt(SupplementaryAlignment::getNumberOfAlignedBases)
                .sum();
    }

    /**
     * @return the readname of this group
     */
    public String getReadName() {
        return original.getReadName();
    }

    /**
     * @return all the unique contigs in this group, sorted
     */
    public List<String> getContigs() {
        return streamInReadOrder()
                .map(SupplementaryAlignment::getContig)
                .filter(Objects::nonNull)
                .sorted(ChromosomeNameComparator.get())
                .distinct()
                .collect(Collectors.toList());
    }


    /**
     * @return the total number of reference bases aligned to this collection of reads
     * this includes deletions but not insertions or clips
     */
    public int getLengthOnReference() {
        return streamInPositionOrder()
                .mapToInt(SupplementaryAlignment::getLengthOnReference)
                .sum();
    }


    public SupplementaryAlignment getNextInRead(SupplementaryAlignment alignment) {
        final int i = readOrder.indexOf(alignment);
        return readOrder.size() > i + 1 ? readOrder.get(i + 1) : null;
    }

    public SupplementaryAlignment getPreviousInRead(SupplementaryAlignment alignment) {
        final int i = readOrder.indexOf(alignment);
        return i > 0 ? readOrder.get(i - 1) : null;
    }

    public SupplementaryAlignment getNextPosition(SupplementaryAlignment alignment) {
        final int i = positionOrder.indexOf(alignment);
        return positionOrder.size() > i + 1 ? readOrder.get(i + 1) : null;
    }

    public SupplementaryAlignment getPreviousPosition(SupplementaryAlignment alignment) {
        final int i = readOrder.indexOf(alignment);
        return i > 0 ? readOrder.get(i - 1) : null;
    }

    /**
     * {@link SAMAlignment} and {@link SupplementaryAlignment} don't share any useful interface. We use
     * and adapter for the read which owns the SA tag (which isn't included in the tag itself) in order to
     * include it seamlessly.
     *
     * @return the adapter wrapping the original read which contains the SA tag this group was built from
     */
    public SupplementaryAlignment getPrimaryAlignment() {
        return primaryRead;
    }

    /**
     * @return the alignment this group was created from
     */
    public Alignment unwrap() {
        return original;
    }

    public Iterator<SupplementaryAlignment> iterateInReadOrder() {
        return readOrder.iterator();
    }

    public Iterator<SupplementaryAlignment> iterateInPositionOrder() {
        return positionOrder.iterator();
    }

    public Stream<SupplementaryAlignment> streamInReadOrder() {
        return readOrder.stream();
    }

    public Stream<SupplementaryAlignment> streamInPositionOrder() {
        return positionOrder.stream();
    }

    public int size() {
        return readOrder.size();
    }

    private static class Adapter extends SupplementaryAlignment {
        public Adapter(Alignment a) {
            super(a.getChr(), a.getStart(), a.getReadStrand(), a.getCigar(), a.getMappingQuality(), a.getAttribute(SAMTag.NM.name()) == null ? 0 : (int) a.getAttribute(SAMTag.NM.name()));
        }
    }
}
