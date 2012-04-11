package org.broad.igv.sam;

import java.util.Arrays;

/**
 * Represents a flow signals context in an alignment block.  Added to support IonTorrent alignments.
 *
 * @author Nils Homer
 * @date 4/11/12
 */
public class FlowSignalContext {
    public short[][][]signals = null;

    public FlowSignalContext(short[][][] signals) {
        this.signals = signals;
    }
}
