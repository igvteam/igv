/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.sam;

import java.util.Arrays;

/**
 * Builds a flow signals context in an alignment block.  Added to support IonTorrent alignments.
 *
 * @author Nils Homer
 * @date 4/11/12
 * Modified by Chantal Roth, 6/21/2012
 */
public class FlowSignalContextBuilder {

    private short[] flowSignals = null;
    private String flowOrder = null;
    private int flowSignalsIndex = -1;
    private int flowOrderIndex = -1;
    private int prevFlowSignalsStart = -1;
    private int prevFlowSignalsEnd = -1;
    private int flowOrderStart = -1;
    private boolean readNegativeStrandFlag;
    private boolean[] incorporations = null; // required for the reverse strand

    
    public static final int PREV = 0;
    public static final int CURR = 1;
    public static final int NEXT = 2;
    
    public FlowSignalContextBuilder(short[] flowSignals, String flowOrder, int flowOrderStart, byte[] readBases, int fromIdx, boolean readNegativeStrandFlag) {
        if (null == flowSignals || null == flowOrder || flowOrderStart < 0) {
            return;
        }

        this.flowSignals = flowSignals;
        this.flowOrder = flowOrder;
        this.flowOrderIndex = this.flowOrderStart = flowOrderStart;
        this.flowSignalsIndex = 0; // NB: the key sequence/barcode sequence should have been removed for the signals already
        this.readNegativeStrandFlag = readNegativeStrandFlag;

        // init
        if (this.readNegativeStrandFlag) {
            int i;
            this.incorporations = new boolean[this.flowSignals.length];
            // go to the end of the signals to find the first sequenced base
            for (i=readBases.length-1;0<=i;i--) {
                while (this.flowOrder.charAt(this.flowOrderIndex) != SAMAlignment.NT2COMP[readBases[i]]) {
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
        /*
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
                                this.flowOrder.charAt(this.flowOrderIndex) != PicardAlignment.NT2COMP[readBases[j]]) { // NB: malicious input can cause infinite loops here
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
        */
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
        
        int[] flowOrderIndices = new int[nBases];
       
        while (0 <= this.flowSignalsIndex && this.flowSignalsIndex < this.flowSignals.length && i < fromIdx + nBases) {
            short s = this.flowSignals[this.flowSignalsIndex];
            char f = this.flowOrder.charAt((this.flowSignalsIndex + this.flowOrderStart) % this.flowOrder.length());
            flowOrderIndices[idx] = flowSignalsIndex+flowOrderStart;
            int nextFlowSignalsStart = -1, nextFlowSignalsEnd = -1;
            int basepos = i + 1;
            if (basepos < readBases.length) {
                if (this.readNegativeStrandFlag) {
                    nextFlowSignalsEnd = this.flowSignalsIndex - 1;
                    // NB: loop condition is not symmetric to the forward, as we must respect the directionality of sequencing.
                    // For example, if our flow order is TACAG, and our read bases are TAG, then the flow signal vector is 
                    // approximately 100,100,0,0,100.  Since we move in the reverse direction with respect to the flow signal 
                    // vector we must pre-compute where the flows incorporations are expected to occur, instead of just looking 
                    // for the next flow that matches our next read base (we would place the A incorporation flow in the fourth flow,
                    // which is wrong).
                    while (!this.incorporations[this.flowSignalsIndex] ||
                            this.flowOrder.charAt(this.flowOrderIndex) != SAMAlignment.NT2COMP[readBases[basepos]]) { // NB: malicious input can cause infinite loops here
                        this.flowOrderIndex--;
                        this.flowSignalsIndex--;
                        if (this.flowOrderIndex < 0) {
                            this.flowOrderIndex = this.flowOrder.length() - 1;
                        }
                            }
                    nextFlowSignalsStart = this.flowSignalsIndex + 1;
                } else {
                    nextFlowSignalsStart = this.flowSignalsIndex + 1;
                    while (this.flowOrder.charAt(this.flowOrderIndex) != readBases[basepos]) { // NB: malicious input can cause infinite loops here
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
                blockFlowSignals[idx][PREV] = new short[this.prevFlowSignalsEnd - this.prevFlowSignalsStart + 1];
                blockFlowOrder[idx][PREV] = new char[this.prevFlowSignalsEnd - this.prevFlowSignalsStart + 1];
                if (this.readNegativeStrandFlag) {
                    for (int flowpos = this.prevFlowSignalsEnd; this.prevFlowSignalsStart <= flowpos; flowpos--) {
                        blockFlowSignals[idx][PREV][this.prevFlowSignalsEnd - flowpos] = this.flowSignals[flowpos];
                        blockFlowOrder[idx][PREV][this.prevFlowSignalsEnd - flowpos] = this.flowOrder.charAt((flowpos + this.flowOrderStart) % this.flowOrder.length());
                    }
                } else {
                    for (int flowpos = this.prevFlowSignalsStart; flowpos <= this.prevFlowSignalsEnd; flowpos++) {
                        blockFlowSignals[idx][PREV][flowpos-this.prevFlowSignalsStart] = this.flowSignals[flowpos];
                        blockFlowOrder[idx][PREV][flowpos-this.prevFlowSignalsStart] = this.flowOrder.charAt((flowpos + this.flowOrderStart) % this.flowOrder.length());
                    }
                }
            } else {
                blockFlowSignals[idx][PREV] = null;
                blockFlowOrder[idx][PREV] = null;
            }
            // current context
            blockFlowSignals[idx][CURR] = new short[1];
            blockFlowOrder[idx][CURR] = new char[1];
            blockFlowSignals[idx][CURR][0] = s;
            blockFlowOrder[idx][CURR][0] = f;
            // next context
            if (0 <= nextFlowSignalsStart && nextFlowSignalsStart <= nextFlowSignalsEnd && nextFlowSignalsEnd < this.flowSignals.length) {
                blockFlowSignals[idx][NEXT] = new short[nextFlowSignalsEnd - nextFlowSignalsStart + 1];
                blockFlowOrder[idx][NEXT] = new char[nextFlowSignalsEnd - nextFlowSignalsStart + 1];
                if (this.readNegativeStrandFlag) {
                    for (int flowpos = nextFlowSignalsEnd; nextFlowSignalsStart <= flowpos; flowpos--) {
                        blockFlowSignals[idx][NEXT][nextFlowSignalsEnd - flowpos] = this.flowSignals[flowpos];
                        blockFlowOrder[idx][NEXT][nextFlowSignalsEnd - flowpos] = this.flowOrder.charAt((flowpos + this.flowOrderStart) % this.flowOrder.length());
                    }
                } else {
                    for (int flowpos = nextFlowSignalsStart; flowpos <= nextFlowSignalsEnd; flowpos++) {
                        blockFlowSignals[idx][NEXT][flowpos-nextFlowSignalsStart] = this.flowSignals[flowpos];
                        blockFlowOrder[idx][NEXT][flowpos-nextFlowSignalsStart] = this.flowOrder.charAt((flowpos + this.flowOrderStart) % this.flowOrder.length());
                    }
                }
            } else {
                blockFlowSignals[idx][NEXT] = null;
                blockFlowOrder[idx][NEXT] = null;
            }
            // update for the next iteration
            this.prevFlowSignalsStart = nextFlowSignalsStart;
            this.prevFlowSignalsEnd = nextFlowSignalsEnd;
            i++; // next base
            idx++; // next base
        }

        return new FlowSignalContext(blockFlowSignals, blockFlowOrder, flowOrderIndices );
    }
}
