package org.igv.sam;

import htsjdk.samtools.util.Locatable;
import org.igv.feature.genome.ChromosomeNameComparator;

import java.util.Comparator;
import java.util.Objects;
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
    }, BASE {
        @Override
        Comparator<Alignment> getAlignmentComparator(final int center, final String tag, final byte referenceBase) {

            final Comparator<Alignment> insertionComparator = Comparator.comparing((Alignment alignment) -> {
                String insertionBases = "";
                AlignmentBlock leftInsertion = alignment.getInsertionAt(center + 1); //todo figure out what's going on with the +1 here..
                if (leftInsertion != null) {
                    insertionBases += leftInsertion.getBases().getString();
                }
                AlignmentBlock rightInsertion = alignment.getInsertionAt(center);
                if (rightInsertion != null) {
                    insertionBases += rightInsertion.getBases().getString();
                }
                return insertionBases;
            }).reversed();

            final int refBase = referenceBase >= 97 ? referenceBase - 32 : referenceBase;

            final ToIntFunction<Alignment> baseCompare = (Alignment a) -> {
                int base = a.getBase(center);
                if (base >= 97) base -= 32; //normalize case
                if (base == refBase) base = Integer.MAX_VALUE - 3;
                if (base == (int) 'N') base = Integer.MAX_VALUE - 2;
                if (base == 0) base = 97;  // > any letter, causes deletions to be placed after snps
                return base;
            };

            final Comparator<Alignment> deletionComparator = Comparator.comparing(
                    (Alignment a) -> a.getDeletionAt(center),
                    Comparator.nullsLast(Comparator.comparing(Gap::getnBases).thenComparing(Gap::getStart)
                    ));

            final ToIntFunction<Alignment> baseQualityCompare = ((Alignment a) -> -a.getPhred(center));

            return Comparator.comparingInt(baseCompare)
                    .thenComparing(deletionComparator)
                    .thenComparing(insertionComparator)
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
                    .thenComparing(a -> a.getMate() != null && Objects.equals(a.getMate().getChr(), a.getChr()))
                    .thenComparing(nullSafeComparator(a -> a.getMate() == null ? null : a.getMate().getChr()));

        }
    }, TAG {
        @Override
        Comparator<Alignment> getAlignmentComparator(final int center, final String tag, final byte referenceBase) {
            return Comparator.comparing((Alignment a) -> a.getAttribute(tag),
                    Comparator.nullsLast(Comparator.comparing(Object::hashCode)));
            //todo It would be nice to sort by something smarter than hash code but the possibility of mixed tag types makes that more complicated
        }
    },
    LEFT_CLIP {
        @Override
        Comparator<Alignment> getAlignmentComparator(final int center, final String tag, final byte referenceBase) {
            return Comparator.comparing((Alignment a) -> a.getClippingCounts().isLeftClipped()).reversed()
                    .thenComparing(Alignment::getAlignmentStart);
        }
    },
    RIGHT_CLIP {
        @Override
        Comparator<Alignment> getAlignmentComparator(final int center, final String tag, final byte referenceBase) {
            return Comparator.comparing((Alignment a) -> a.getClippingCounts().isRightClipped()).reversed()
                    .thenComparing(Alignment::getAlignmentEnd);
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
            return Comparator.comparingInt(Alignment::getClusterDistance);
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
            return Comparator.comparingInt(a -> a.getAlignmentStart() - a.getAlignmentEnd());
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
        if (invertSort) alignmentComparator = alignmentComparator.reversed();
        return Comparator.comparing((Row row) -> AlignmentInterval.getFeatureContaining(row.getAlignments(), center),
                Comparator.nullsLast(alignmentComparator));

    }

    // private method to save a few lines of code in each comparator since they all would have to start with getting the
    // center alignment anyway
    abstract Comparator<Alignment> getAlignmentComparator(final int center, final String tag, final byte referenceBase);

    /**
     * Custom valueOf method with backward compatibility.
     * Supports "NUCLEOTIDE" as an alias for "BASE" for backward compatibility.
     *
     * @param name the string name of the enum constant
     * @return the enum constant with the specified name
     * @throws IllegalArgumentException if no constant with the specified name is found
     */
    public static SortOption fromString(String name) {
        if (name == null) {
            return SortOption.NONE;
        }
        if ("NUCLEOTIDE".equals(name)) {
            return BASE;
        }
        return SortOption.valueOf(name);
    }

    public static final Comparator<Locatable> POSITION_COMPARATOR = Comparator.nullsFirst(Comparator.comparing(Locatable::getContig, ChromosomeNameComparator.get()))
            .thenComparing(Locatable::getStart)
            .thenComparing(Locatable::getEnd);
}
