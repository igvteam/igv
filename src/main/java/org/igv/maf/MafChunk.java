package org.igv.maf;

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
