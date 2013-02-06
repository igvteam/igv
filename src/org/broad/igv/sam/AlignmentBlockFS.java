/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

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

    protected AlignmentBlockFS(int start, byte[] bases, byte[] qualities, FlowSignalContext fContext) {
        super(start, bases, qualities);
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
