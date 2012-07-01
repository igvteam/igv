package org.broad.igv.sam;

/**
 * Represents an alignment block which contains flow signals.  Added to support IonTorrent alignments.
 *
 * @author Jim Robinson
 * @date 3/19/12
 * Modified by Chantal Roth, 6/21/2012
 */
public class AlignmentBlockFS extends AlignmentBlock {

    private FlowSignalContext fContext = null;

    protected AlignmentBlockFS(int start, byte[] bases, byte[] qualities, FlowSignalContext fContext, Alignment baseAlignment) {
        super(start, bases, qualities, baseAlignment);
        if (fContext != null && fContext.getNrSignals() == bases.length) {
            this.fContext = fContext;
        }
    }
    public FlowSignalContext getFlowSignalContext() {
        return fContext;
    }

    @Override
    public FlowSignalSubContext getFlowSignalSubContext(int offset) {
        return new FlowSignalSubContext(this.fContext.getSignalForOffset(offset), this.fContext.getBasesForOffset(offset), this.fContext.getFlowOrderIndexForOffset(offset));
    }


    @Override
    public boolean hasFlowSignals() {
        return (null != this.fContext);
    }
}
