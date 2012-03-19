package org.broad.igv.sam;

/**
 * Represents an alignment block which contains flow signals.  Added to suppor IonTorrent alignments.
 *
 * @author Jim Robinson
 * @date 3/19/12
 */
public class AlignmentBlockFS extends AlignmentBlock {

    public short[][][] flowSignals = null;

    protected AlignmentBlockFS(int start, byte[] bases, byte[] qualities, short[][][] flowSignals, Alignment baseAlignment) {
        super(start, bases, qualities, baseAlignment);
        if (flowSignals != null && flowSignals.length == bases.length) {
            this.flowSignals = flowSignals;
        }
    }

    public short[][] getFlowSignalContext(int offset) {
        return this.flowSignals[offset];
    }


    public boolean hasFlowSignals() {
        return (null != this.flowSignals);
    }
}
