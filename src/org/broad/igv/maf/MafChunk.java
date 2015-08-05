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

package org.broad.igv.maf;

import java.util.*;

/**
 * @author jrobinso
 *         Date: 2/19/13
 *         Time: 1:09 AM
 */
public class MafChunk {

    String chr;
    int start;
    int end;
    LinkedHashMap<String, MultipleAlignmentBlock> blocks;


    public MafChunk(String chr, int start, int end, List<MultipleAlignmentBlock> blockList) {
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.blocks = new LinkedHashMap<String, MultipleAlignmentBlock>();
        for (MultipleAlignmentBlock b : blockList) {
            blocks.put(b.getKey(), b);
        }
    }


    public boolean overlaps(String chr, int start, int end) {

        if(!chr.equals(this.chr)) return false;
        return end >= this.start && start <= this.end;
    }

    public boolean contains(String chr, int start, int end) {
        if(!chr.equals(this.chr)) return false;
        return start >= this.start && end <= this.end;

    }


    public void addBlocks(int start, int end, List<MultipleAlignmentBlock> blockList) {

        // Can only add blocks if interval is contiguous with current extent.  Otherwise start over.
        if (start > this.end || end < this.start) {
            blocks.clear();
            for (MultipleAlignmentBlock b : blockList) {
                blocks.put(b.getKey(), b);
            }
        } else if (start < this.start) {
            // Add to left
            LinkedHashMap<String, MultipleAlignmentBlock> tmp = this.blocks;
            this.blocks = new LinkedHashMap<String, MultipleAlignmentBlock>();
            for (MultipleAlignmentBlock b : blockList) {
                String key = b.getKey();
                if(blocks.containsKey(key)) {
                    // Done
                    break;
                }
                blocks.put(key, b);
            }
            this.blocks.putAll(tmp);
            this.start = start;
        } else {
            // Add to right
            for (MultipleAlignmentBlock b : blockList) {
                String key = b.getKey();
                if(blocks.containsKey(key)) {
                    continue;
                }
                blocks.put(b.getKey(), b);
            }
            this.end = end;
        }
    }

    public void trimTo(int start, int end) {

        LinkedHashMap<String, MultipleAlignmentBlock> tmp = blocks;
        blocks = new LinkedHashMap<String, MultipleAlignmentBlock>();
        for(Map.Entry<String, MultipleAlignmentBlock> entry : tmp.entrySet()) {
            String key = entry.getKey();
            MultipleAlignmentBlock block = entry.getValue();
            if(block.getEnd() < start) {
                continue;
            }
            else if(block.getStart() > end) {
                break;
            }
            else {
                blocks.put(key, block);
            }
        }

    }
}
