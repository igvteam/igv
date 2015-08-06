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

import org.broad.igv.AbstractHeadlessTest;
import org.junit.Ignore;

/**
 * User: jacob
 * Date: 2013-Feb-05
 */
@Ignore("No runnable tests")
public class AlignmentBlockTest extends AbstractHeadlessTest{

//    @Test
//    public void testCreateMismatchBlocks() throws Exception{
//        AlignmentDataManager manager = AlignmentDataManagerTest.getManager171();
//
//        final String chr = "chr1";
//        final int queryStart = 151666494;
//        final int halfwidth = 1000;
//        final int queryEnd = queryStart + 2 * halfwidth;
//
//        AlignmentInterval interval = manager.loadInterval(chr, queryStart, queryEnd, new AlignmentTrack.RenderOptions());
//
//        Iterator<Alignment> alignmentIterator = interval.getAlignmentIterator();
//        while(alignmentIterator.hasNext()){
//            Alignment alignment = alignmentIterator.next();
//            AlignmentBlock[] blocks = alignment.getAlignmentBlocks();
//            for(AlignmentBlock block: blocks){
//
//                int start = block.getStart();
//                int end = block.getEnd();
//                byte[] refSeq = genome.getSequence(chr, start, end);
//
//                AlignmentBlock.MismatchBlock[] mismatchBlocks = AlignmentBlock.createMismatchBlocks(start, refSeq, block.getBases());
//                AlignmentBlock.MismatchBlock nextBlock = null;
//                int blockInd = 0;
//                int nextBlockStart = Integer.MAX_VALUE;
//                int nextBlockEnd = Integer.MIN_VALUE;
//                int mmIdx = 0;
//
//                if(mismatchBlocks.length > 0){
//                    nextBlock = mismatchBlocks[blockInd++];
//                    nextBlockStart = nextBlock.start;
//                    nextBlockEnd = nextBlockStart + nextBlock.bases.length;
//                }
//
//                for(int loc = start; loc < end; loc++){
//
//                    int idx = loc - start;
//                    byte base = block.getBase(idx);
//
//                    if(loc >= nextBlockStart && loc < nextBlockEnd){
//                        assertEquals(base, nextBlock.bases[mmIdx++]);
//                    }else if(loc == nextBlockEnd){
//                        if(blockInd < mismatchBlocks.length){
//                            nextBlock = mismatchBlocks[blockInd++];
//                            nextBlockStart = nextBlock.start;
//                            nextBlockEnd = nextBlockStart + nextBlock.bases.length;
//                            mmIdx = 0;
//                        }else{
//                            nextBlock = null;
//                            nextBlockStart = Integer.MAX_VALUE;
//                            nextBlockEnd = Integer.MIN_VALUE;
//                        }
//                    }else{
//                        String msg = String.format("Expected match between read and reference at %d, but was mismatch: %s, %s. Block: %s:%s-%s", idx, (char) refSeq[idx], (char) base, chr, start+1, end);
//                        assertTrue(msg, AlignmentUtils.compareBases(refSeq[idx], base));
//                    }
//                }
//            }
//        }
//
//    }

//    @Ignore
//    @Test
//    public void testMemoryFootprintMismatchBlock() throws Exception{
//        System.gc();
//        int numMis = 1000000;
//        long before = RuntimeUtils.getAvailableMemory();
//
//        AlignmentBlock.MismatchBlock[] blocks = new AlignmentBlock.MismatchBlock[numMis];
//        for(int ii=0; ii < numMis; ii++){
//            blocks[ii] = new AlignmentBlock.MismatchBlock(0, "thisismaoeuystring".getBytes());
//        }
//        System.gc();
//        long after = RuntimeUtils.getAvailableMemory();
//        float measMem = (float) RuntimeUtils.getObjectSizeRecursive(blocks, new HashSet<Object>());
//
//        float diff = (float) (before - after);
//        System.out.println(String.format("Difference before/after: %f kb. %d objects. %f bytes each", (before - after) / 1000.0, numMis, diff/numMis));
//        System.out.println(String.format("Measured memory: %f kb. %d objects. %f bytes each", measMem / 1000.0, numMis, measMem / numMis));
//    }


}
