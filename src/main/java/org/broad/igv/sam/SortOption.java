package org.broad.igv.sam;

import java.util.Comparator;
import java.util.function.Function;
import java.util.function.ToIntFunction;

public enum SortOption {
    START {
        @Override
        Comparator<Alignment> getAlignmentComparator(final int center, final String tag, final byte referenceBase) {
            return Comparator.comparingInt(Alignment::getAlignmentStart);
        }
    }, STRAND {
        @Override
        Comparator<Alignment> getAlignmentComparator(final int center, final String tag, final byte referenceBase) {
            return nullSafeComparator((Alignment a) -> a instanceof LinkedAlignment
                    ? ((LinkedAlignment) a).getStrandAtPosition(center)
                    : a.getReadStrand());
        }
    }, NUCLEOTIDE {
        @Override
        Comparator<Alignment> getAlignmentComparator(final int center, final String tag, final byte referenceBase) {

            final Comparator<Alignment> insertionComparator = Comparator.comparing((Alignment alignment) -> {
                AlignmentBlock insertion = alignment.getInsertionAt(center + 1); //todo figure out what's going on with the +1 here..
                if(insertion == null) insertion = alignment.getInsertionAt(center);
                return insertion != null ? insertion.getBases().getString() : "";
            }).reversed();

            final int refBase = referenceBase >= 97 ? referenceBase - 32 : referenceBase;

            final ToIntFunction<Alignment> baseCompare = (Alignment a) -> {
                int base = a.getBase(center);
                if (base >= 97) base -= 32; //normalize case
                if (base == refBase) base = Integer.MAX_VALUE - 3;
                if (base == (int)'N') base = Integer.MAX_VALUE - 2;
                if (base == 0) base = Integer.MAX_VALUE;
                return base;
            };

            final Comparator<Alignment> deletionComparator = Comparator.comparing(
                        (Alignment a) -> a.getDeletionAt(center),
                        Comparator.nullsLast(Comparator.comparing(Gap::getnBases).thenComparing(Gap::getStart)
                    ));

            final ToIntFunction<Alignment> baseQualityCompare = ((Alignment a) -> -a.getPhred(center));

            return insertionComparator
                    .thenComparing(deletionComparator)
                    .thenComparingInt(baseCompare)
                    .thenComparingInt(baseQualityCompare);
        }
    }, QUALITY {
        @Override
        Comparator<Alignment> getAlignmentComparator(final int center, final String tag, final byte referenceBase) {
            return Comparator.comparingInt(Alignment::getMappingQuality).reversed();
        }
    }, SAMPLE {
        @Override
        Comparator<Alignment> getAlignmentComparator(final int center, final String tag, final byte referenceBase) {
            return nullSafeComparator(Alignment::getSample);
        }
    }, READ_GROUP {
        @Override
        Comparator<Alignment> getAlignmentComparator(final int center, final String tag, final byte referenceBase) {
            return nullSafeComparator(Alignment::getReadGroup);
        }
    }, INSERT_SIZE {
        @Override
        Comparator<Alignment> getAlignmentComparator(final int center, final String tag, final byte referenceBase) {
            return Comparator.comparingInt((Alignment a) -> Math.abs(a.getInferredInsertSize())).reversed();
        }
    }, FIRST_OF_PAIR_STRAND {
        @Override
        Comparator<Alignment> getAlignmentComparator(final int center, final String tag, final byte referenceBase) {
            return nullSafeComparator(Alignment::getFirstOfPairStrand);
        }
    }, MATE_CHR {
        @Override
        Comparator<Alignment> getAlignmentComparator(final int center, final String tag, final byte referenceBase) {
                return Comparator.comparing((Alignment a) -> a.getMate() == null)
                    .thenComparing(a -> a.getMate().getChr().equals(a.getChr()))
                    .thenComparing(a -> a.getMate().getChr());
        }
    }, TAG {
        @Override
        Comparator<Alignment> getAlignmentComparator(final int center, final String tag, final byte referenceBase) {
            return Comparator.comparing((Alignment a) -> a.getAttribute(tag),
                    Comparator.nullsLast(Comparator.comparing(Object::hashCode)));
            //todo It would be nice to sort by something smarter than hash code but the possibility of mixed tag types makes that more complicated
        }
    },
    SUPPLEMENTARY {
        @Override
        Comparator<Alignment> getAlignmentComparator(final int center, final String tag, final byte referenceBase) {
            return Comparator.comparing(Alignment::isSupplementary).reversed();
        }
    }, NONE {
        @Override
        Comparator<Alignment> getAlignmentComparator(final int center, final String tag, final byte referenceBase) {
            return (a1, a2) -> 0;
        }
    }, HAPLOTYPE {
        @Override
        Comparator<Alignment> getAlignmentComparator(final int center, final String tag, final byte referenceBase) {
            return Comparator.comparingInt(Alignment::getHapDistance);
        }
    }, READ_ORDER {
        @Override
        Comparator<Alignment> getAlignmentComparator(final int center, final String tag, final byte referenceBase) {
            return Comparator.comparing(Alignment::isPaired)
                    .thenComparing(Alignment::isFirstOfPair)
                    .thenComparing(Alignment::isSecondOfPair).reversed();
        }
    }, READ_NAME {
        @Override
        Comparator<Alignment> getAlignmentComparator(final int center, final String tag, final byte referenceBase) {
            return Comparator.comparing(Alignment::getReadName);
        }
    }, ALIGNED_READ_LENGTH {
        @Override
        Comparator<Alignment> getAlignmentComparator(final int center, final String tag, final byte referenceBase) {
            return Comparator.comparingInt( a -> a.getAlignmentStart() - a.getAlignmentEnd());
        }
    };

    // Slightly shortened way of writing this common use case
    protected final <T, U extends Comparable<? super U>> Comparator<T> nullSafeComparator(final Function<? super T, ? extends U> extractor) {
        return Comparator.comparing(extractor, Comparator.nullsLast(Comparator.naturalOrder()));
    }

    /**
     * get a Comparator that will perform the relevant sort. 
     */
    public Comparator<Row> getComparator(final int center, final byte reference, final String tag, boolean invertSort) {
        Comparator<Alignment> alignmentComparator = getAlignmentComparator(center, tag, reference);
        if(invertSort) alignmentComparator = alignmentComparator.reversed();
        return Comparator.comparing((Row row) -> AlignmentInterval.getFeatureContaining(row.getAlignments(), center),
                Comparator.nullsLast(alignmentComparator));

    }

    // private method to save a few lines of code in each comparator since they all would have to start with getting the
    // center alignment anyway
    abstract Comparator<Alignment> getAlignmentComparator(final int center, final String tag, final byte referenceBase);

}
