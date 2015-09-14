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
