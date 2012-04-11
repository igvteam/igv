package org.broad.igv.sam;

import java.util.Arrays;

/**
 * Builds a flow signals context in an alignment block.  Added to support IonTorrent alignments.
 *
 * @author Nils Homer
 * @date 4/11/12
 */
public class FlowSignalContextBuilder {

    private short[] flowSignals = null;
    private String flowOrder = null;
    private int flowSignalsIndex = -1;
    private int flowOrderIndex = -1;
    private int prevFlowSignalsStart = -1;
    private int prevFlowSignalsEnd = -1;
    private boolean readNegativeStrandFlag;
    private boolean[] incorporations = null; // required for the reverse strand

    public FlowSignalContextBuilder(short[] flowSignals, String flowOrder, int flowOrderStart, byte[] readBases, int fromIdx, boolean readNegativeStrandFlag) {
        if (null == flowSignals || null == flowOrder || flowOrderStart < 0) {
            return;
        }

        this.flowSignals = flowSignals;
        this.flowOrder = flowOrder;
        this.flowOrderIndex = flowOrderStart;
        this.flowSignalsIndex = 0; // NB: the key sequence/barcode sequence should have been remove for the signals already
        this.readNegativeStrandFlag = readNegativeStrandFlag;

        // init
        if (this.readNegativeStrandFlag) {
            int i;
            this.incorporations = new boolean[this.flowSignals.length];
            // go to the end of the signals
            for (i=readBases.length-1;0<=i;i--) {
                while (this.flowOrder.charAt(this.flowOrderIndex) != SamAlignment.NT2COMP[readBases[i]]) {
                    this.flowOrderIndex++;
                    this.flowSignalsIndex++;
                    this.incorporations[this.flowSignalsIndex] = false;
                    if (this.flowOrder.length() <= this.flowOrderIndex) {
                        this.flowOrderIndex = 0;
                    }
                }
                this.incorporations[this.flowSignalsIndex] = true;
            }
            this.prevFlowSignalsStart = this.flowSignalsIndex + 1;
            this.prevFlowSignalsEnd = this.flowSignals.length - 1;
        } else {
            this.prevFlowSignalsStart = this.prevFlowSignalsEnd = 0;
            while (this.flowOrder.charAt(this.flowOrderIndex) != readBases[0]) {
                this.flowOrderIndex++;
                this.flowSignalsIndex++;
                if (this.flowOrder.length() <= this.flowOrderIndex) {
                    this.flowOrderIndex = 0;
                }
            }
            this.prevFlowSignalsEnd = this.flowSignalsIndex - 1;
        }
        if (0 < fromIdx) { // skip over leading bases (ex. soft clipped bases)
            int i = 0;
            while (0 <= this.flowSignalsIndex && this.flowSignalsIndex < this.flowSignals.length && i < fromIdx) {
                short s = this.flowSignals[this.flowSignalsIndex];
                int nextFlowSignalsStart = -1, nextFlowSignalsEnd = -1;
                int j = i + 1;
                if (j < readBases.length) {
                    if (this.readNegativeStrandFlag) {
                        nextFlowSignalsEnd = this.flowSignalsIndex - 1;
                        // NB: loop condition is not symmetric to the forward, as we must respect the directionality of sequencing.
                        // For example, if our flow order is TACAG, and our read bases are TAG, then the flow signal vector is 
                        // approximately 100,100,0,0,100.  Since we move in the reverse direction with respect to the flow signal 
                        // vector we must pre-compute where the flows incorporations are expected to occur, instead of just looking 
                        // for the next flow that matches our next read base (we would place the A incorporation flow in the fourth flow,
                        // which is wrong).
                        while (!this.incorporations[this.flowSignalsIndex] ||
                                this.flowOrder.charAt(this.flowOrderIndex) != SamAlignment.NT2COMP[readBases[j]]) { // NB: malicious input can cause infinite loops here
                            this.flowOrderIndex--;
                            this.flowSignalsIndex--;
                            if (this.flowOrderIndex < 0) {
                                this.flowOrderIndex = this.flowOrder.length() - 1;
                            }
                                }
                        nextFlowSignalsStart = this.flowSignalsIndex + 1;
                    } else {
                        nextFlowSignalsStart = this.flowSignalsIndex + 1;
                        while (this.flowOrder.charAt(this.flowOrderIndex) != readBases[j]) { // NB: malicious input can cause infinite loops here
                            this.flowOrderIndex++;
                            this.flowSignalsIndex++;
                            if (this.flowOrder.length() <= this.flowOrderIndex) {
                                this.flowOrderIndex = 0;
                            }
                        }
                        nextFlowSignalsEnd = this.flowSignalsIndex - 1;
                    }
                }
                // update for the next iteration
                this.prevFlowSignalsStart = nextFlowSignalsStart;
                this.prevFlowSignalsEnd = nextFlowSignalsEnd;
                i++; // next base
            }
        }
    }

    // TODO:
    // - support IUPAC bases
    // - support lower/upper cases (is this necessary)?
    public FlowSignalContext getFlowSignalContext(byte[] readBases, int fromIdx, int nBases) {
        int i, idx;
        short[][][] blockFlowSignals = null;
        char[][][] blockFlowOrder = null; 

        if (null == this.flowSignals) {
            return null;
        }

        blockFlowSignals = new short[nBases][][];
        blockFlowOrder = new char[nBases][][];
        //Default value
        Arrays.fill(blockFlowSignals, null);
        Arrays.fill(blockFlowOrder, null);

        // NB: should be at the first base of a HP
        // Go through the bases
        i = fromIdx;
        idx = 0;
        while (0 <= this.flowSignalsIndex && this.flowSignalsIndex < this.flowSignals.length && i < fromIdx + nBases) {
            short s = this.flowSignals[this.flowSignalsIndex];
            char f = this.flowOrder.charAt(this.flowSignalsIndex % this.flowOrder.length());
            int nextFlowSignalsStart = -1, nextFlowSignalsEnd = -1;
            int j = i + 1;
            if (j < readBases.length) {
                if (this.readNegativeStrandFlag) {
                    nextFlowSignalsEnd = this.flowSignalsIndex - 1;
                    // NB: loop condition is not symmetric to the forward, as we must respect the directionality of sequencing.
                    // For example, if our flow order is TACAG, and our read bases are TAG, then the flow signal vector is 
                    // approximately 100,100,0,0,100.  Since we move in the reverse direction with respect to the flow signal 
                    // vector we must pre-compute where the flows incorporations are expected to occur, instead of just looking 
                    // for the next flow that matches our next read base (we would place the A incorporation flow in the fourth flow,
                    // which is wrong).
                    while (!this.incorporations[this.flowSignalsIndex] ||
                            this.flowOrder.charAt(this.flowOrderIndex) != SamAlignment.NT2COMP[readBases[j]]) { // NB: malicious input can cause infinite loops here
                        this.flowOrderIndex--;
                        this.flowSignalsIndex--;
                        if (this.flowOrderIndex < 0) {
                            this.flowOrderIndex = this.flowOrder.length() - 1;
                        }
                            }
                    nextFlowSignalsStart = this.flowSignalsIndex + 1;
                } else {
                    nextFlowSignalsStart = this.flowSignalsIndex + 1;
                    while (this.flowOrder.charAt(this.flowOrderIndex) != readBases[j]) { // NB: malicious input can cause infinite loops here
                        this.flowOrderIndex++;
                        this.flowSignalsIndex++;
                        if (this.flowOrder.length() <= this.flowOrderIndex) {
                            this.flowOrderIndex = 0;
                        }
                    }
                    nextFlowSignalsEnd = this.flowSignalsIndex - 1;
                }
            }
            // set-up block
            blockFlowSignals[idx] = new short[3][];
            blockFlowOrder[idx] = new char[3][];
            // this.previous context
            if (0 <= this.prevFlowSignalsStart && this.prevFlowSignalsStart <= this.prevFlowSignalsEnd && this.prevFlowSignalsEnd < this.flowSignals.length) {
                blockFlowSignals[idx][0] = new short[this.prevFlowSignalsEnd - this.prevFlowSignalsStart + 1];
                blockFlowOrder[idx][0] = new char[this.prevFlowSignalsEnd - this.prevFlowSignalsStart + 1];
                if (this.readNegativeStrandFlag) {
                    for (j = this.prevFlowSignalsEnd; this.prevFlowSignalsStart <= j; j--) {
                        blockFlowSignals[idx][0][this.prevFlowSignalsEnd - j] = this.flowSignals[j];
                        blockFlowOrder[idx][0][this.prevFlowSignalsEnd - j] = this.flowOrder.charAt(j % this.flowOrder.length());
                    }
                } else {
                    for (j = this.prevFlowSignalsStart; j <= this.prevFlowSignalsEnd; j++) {
                        blockFlowSignals[idx][0][j-this.prevFlowSignalsStart] = this.flowSignals[j];
                        blockFlowOrder[idx][0][j-this.prevFlowSignalsStart] = this.flowOrder.charAt(j % this.flowOrder.length());
                    }
                }
            } else {
                blockFlowSignals[idx][0] = null;
                blockFlowOrder[idx][0] = null;
            }
            // current context
            blockFlowSignals[idx][1] = new short[1];
            blockFlowOrder[idx][1] = new char[1];
            blockFlowSignals[idx][1][0] = s;
            blockFlowOrder[idx][1][0] = f;
            // next context
            if (0 <= nextFlowSignalsStart && nextFlowSignalsStart <= nextFlowSignalsEnd && nextFlowSignalsEnd < this.flowSignals.length) {
                blockFlowSignals[idx][2] = new short[nextFlowSignalsEnd - nextFlowSignalsStart + 1];
                blockFlowOrder[idx][2] = new char[nextFlowSignalsEnd - nextFlowSignalsStart + 1];
                if (this.readNegativeStrandFlag) {
                    for (j = nextFlowSignalsEnd; nextFlowSignalsStart <= j; j--) {
                        blockFlowSignals[idx][2][nextFlowSignalsEnd - j] = this.flowSignals[j];
                        blockFlowOrder[idx][2][nextFlowSignalsEnd - j] = this.flowOrder.charAt(j % this.flowOrder.length());
                    }
                } else {
                    for (j = nextFlowSignalsStart; j <= nextFlowSignalsEnd; j++) {
                        blockFlowSignals[idx][2][j-nextFlowSignalsStart] = this.flowSignals[j];
                        blockFlowOrder[idx][2][j-nextFlowSignalsStart] = this.flowOrder.charAt(j % this.flowOrder.length());
                    }
                }
            } else {
                blockFlowSignals[idx][2] = null;
                blockFlowOrder[idx][2] = null;
            }
            // update for the next iteration
            this.prevFlowSignalsStart = nextFlowSignalsStart;
            this.prevFlowSignalsEnd = nextFlowSignalsEnd;
            i++; // next base
            idx++; // next base
        }

        return new FlowSignalContext(blockFlowSignals, blockFlowOrder);
    }
}
