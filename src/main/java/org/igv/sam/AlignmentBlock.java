package org.igv.sam;

import org.igv.feature.genome.Genome;

/**
 * Created by jrobinso on 9/22/15.
 */
public interface AlignmentBlock {

    boolean contains(int position);

    int getBasesLength();

    /**
     * Offset into read sequence
     */
    int getBasesOffset();

    int getLength();

    byte getBase(int offset);

    ByteSubarray getBases();

    int getStart();

    default boolean hasQualities() {
        return false;
    }

    byte getQuality(int offset);

    ByteSubarray getQualities();

    int getEnd();

    boolean isSoftClip();

    boolean hasBases();

    void setPixelRange(int s, int e);

    boolean containsPixel(int x) ;

    default int getPadding() {
        return 0;
    }

    char getCigarOperator();

    default int getIndexOnRead() { return 0; }
}
