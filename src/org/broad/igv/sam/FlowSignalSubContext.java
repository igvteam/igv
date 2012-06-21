package org.broad.igv.sam;

/**
 * Represents a flow signals context in an alignment block focused on a given base.  Added to support IonTorrent alignments.
 *
 * @author Nils Homer
 * @date 4/11/12
 * Modified by Chantal Roth, 6/21/2012
 */
public class FlowSignalSubContext {
    private short[][] signals = null;
    private char[][] bases = null;
    private int flowOrderIndex;
    
    public static final int PREV = 0;
    public static final int CURR = 1;
    public static final int NEXT = 2;
    
    public FlowSignalSubContext(short[][] signals, char[][] bases, int flowOrderIndex) {
        this.signals = signals;
        this.bases = bases;
        this.flowOrderIndex = flowOrderIndex;
    }   
    public int getFlowOrderIndex() {
        return flowOrderIndex;
    }
    public short[] getPreviousSignals() {
        return signals[PREV];
    }
    public short[] getCurrentSignals() {
        return signals[CURR];
    }
     public short getCurrentSignal() {
        return signals[CURR][0];
    }
     public short[] getPNextSignals() {
        return signals[NEXT];
    }

    public int getNrSignalTypes() {
        return signals.length;
    }

    public short[] getSignalsOfType(int type) {
        return signals[type];
    }

    public char[] getBasesOfType(int type) {
        return bases[type];
    }

    public short[][] getSignals() {
        return signals;
    }
     public char[][] getBases() {
        return bases;
    }
} 
