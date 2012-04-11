package org.broad.igv.sam;

import java.util.Arrays;

/**
 * Represents an alignment block which contains flow signals.  Added to support IonTorrent alignments.
 *
 * @author Jim Robinson
 * @date 3/19/12
 */
public class AlignmentBlockFS extends AlignmentBlock {

    public FlowSignalContext fContext = null;

    protected AlignmentBlockFS(int start, byte[] bases, byte[] qualities, FlowSignalContext fContext, Alignment baseAlignment) {
        super(start, bases, qualities, baseAlignment);
        if (fContext != null && fContext.signals.length == bases.length) {
            this.fContext = fContext;
        }
    }

    public short[][] getFlowSignalContext(int offset) {
        return this.fContext.signals[offset];
    }


    public boolean hasFlowSignals() {
        return (null != this.fContext);
    }
}
