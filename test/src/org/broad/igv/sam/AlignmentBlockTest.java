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

import org.broad.igv.AbstractHeadlessTest;
import org.junit.Test;

import java.util.Iterator;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2013-Feb-05
 */
public class AlignmentBlockTest extends AbstractHeadlessTest{

    @Test
    public void testCreateMismatchBlocks() throws Exception{
        AlignmentDataManager manager = AlignmentDataManagerTest.getManager171();

        final String chr = "chr1";
        final int queryStart = 151666494;
        final int halfwidth = 1000;
        final int queryEnd = queryStart + 2 * halfwidth;

        AlignmentInterval interval = manager.loadInterval(chr, queryStart, queryEnd, new AlignmentTrack.RenderOptions());

        Iterator<Alignment> alignmentIterator = interval.getAlignmentIterator();
        while(alignmentIterator.hasNext()){
            Alignment alignment = alignmentIterator.next();
            AlignmentBlock[] blocks = alignment.getAlignmentBlocks();
            for(AlignmentBlock block: blocks){

                int start = block.getStart();
                int end = block.getEnd();
                byte[] refSeq = genome.getSequence(chr, start, end);

                AlignmentBlock.MismatchBlock[] mismatchBlocks = AlignmentBlock.createMismatchBlocks(start, refSeq, block.getBases(), block.getQualities());
                AlignmentBlock.MismatchBlock nextBlock = null;
                int blockInd = 0;
                int nextBlockStart = Integer.MAX_VALUE;
                int nextBlockEnd = Integer.MIN_VALUE;
                int mmIdx = 0;

                if(mismatchBlocks.length > 0){
                    nextBlock = mismatchBlocks[blockInd++];
                    nextBlockStart = nextBlock.start;
                    nextBlockEnd = nextBlockStart + nextBlock.bases.length;
                }

                for(int loc = start; loc < end; loc++){

                    int idx = loc - start;
                    byte base = block.getBase(idx);

                    if(loc >= nextBlockStart && loc < nextBlockEnd){
                        assertEquals(base, nextBlock.bases[mmIdx++]);
                    }else if(loc == nextBlockEnd){
                        if(blockInd < mismatchBlocks.length){
                            nextBlock = mismatchBlocks[blockInd++];
                            nextBlockStart = nextBlock.start;
                            nextBlockEnd = nextBlockStart + nextBlock.bases.length;
                            mmIdx = 0;
                        }else{
                            nextBlock = null;
                            nextBlockStart = Integer.MAX_VALUE;
                            nextBlockEnd = Integer.MIN_VALUE;
                        }
                    }else{
                        String msg = String.format("Expected match between read and reference at %d, but was mismatch: %s, %s. Block: %s:%s-%s", idx, (char) refSeq[idx], (char) base, chr, start+1, end);
                        assertTrue(msg, AlignmentUtils.compareBases(refSeq[idx], base));
                    }
                }
            }
        }

    }


}
