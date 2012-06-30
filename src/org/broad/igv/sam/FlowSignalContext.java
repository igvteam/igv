package org.broad.igv.sam;

import java.util.Arrays;

/**
 * Represents a flow signals context in an alignment block.  Added to support IonTorrent alignments.
 *
 * @author Nils Homer
 * @date 4/11/12
 * modified by Chantal Roth
 */
public class FlowSignalContext {
    private short[][][]signals = null;
    private char[][][]bases = null;
    private int[] flowOrderIndices;

    public FlowSignalContext(short[][][] signals, char[][][] bases, int[] flowOrderIndices) {
        this.signals = signals;
        this.bases = bases;
        this.flowOrderIndices = flowOrderIndices;
    }
    public int getNrSignals() {
        return signals.length;
    }
    public int getNrBases() {
        return bases.length;
    }
    public short[][] getSignalForOffset(int offset) {
        return signals[offset];
    }
    
    public char[][] getBasesForOffset(int offset) {
        return bases[offset];
    }
    public int getFlowOrderIndexForOffset(int offset) {
        return flowOrderIndices[offset];
    }
}
