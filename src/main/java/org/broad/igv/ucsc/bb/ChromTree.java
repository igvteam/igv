package org.broad.igv.ucsc.bb;

import htsjdk.samtools.seekablestream.SeekableStream;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.ucsc.BPTree;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

/**
 * Represents the ChromTree of a UCSC bigbed/bigwig file.  The entire tree is walked to produce 2 maps,
 * (1) ID -> chromosome names, and its
 * (2) chromsome name -> ID
 * <p>
 * Both maps are needed by IGV
 * <p>
 * The chromosome tree is a B+ index, but is located continguously in memory in the header section of the file. This
 * makes it feasible to parse the whole tree with data from a single fetch.   In the end the tree is discarded
 * leaving only the mapps.
 */

public class ChromTree {

    private static Logger log = LogManager.getLogger(ChromTree.class);

    private final long startOffset;
    private BPTree bpTree;
    private HashMap<String, Integer> nameToId = new HashMap<>();
    private HashMap<Integer, String> idToName = new HashMap<>();

    public ChromTree(String file, long startOffset) throws IOException {

        if (FileUtils.isRemote(file)) {
            SeekableStream stream = IGVSeekableStreamFactory.getInstance().getBufferedStream(
                    IGVSeekableStreamFactory.getInstance().getStreamFor(file), 64 * 1024);
            this.bpTree = new BPTree(stream, startOffset);
        } else {
            this.bpTree = new BPTree(file, startOffset);
        }
        this.startOffset = startOffset;
    }

    /**
     * Return the chromosome ID for the given name.  This is the internal chromosome ID for the parent BB file
     * only.
     * @param chr
     * @return
     */
    public Integer getIdForName(String chr) {
        if (nameToId.containsKey(chr)) {
            return nameToId.get(chr);
        } else {
            try {
                int[] result = this.bpTree.searchIntInt(chr);
                if (result != null) {
                    int id = result[0];
                    nameToId.put(chr, id);
                    return id;
                } else {
                    return null;
                }
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
    }

    /**
     * Return the chromosome name for the given ID.  This is a potentially expensive operation as it involves
     * walking the tree until the leaf item for the given name is found.  Currently it is used in only 2
     * situations:
     * (1) decoding features from a bigbed search-by-name query
     * (2) decoding bigwig data from the whole genome view
     *
     * @param id
     * @return
     */
    public String getNameForId(int id) {
        if (idToName.containsKey(id)) {
            return idToName.get(id);
        } else {
            String name = this.searchForName(id);
            if (name != null) {
                idToName.put(id, name);
                return name;
            }
        }
        return null;
    }

    public long getItemCount() {
        return this.bpTree.itemCount;
    }

    /**
     * Perform a reverse search by traversing the tree starting at the given offset.  This is potentially expensive.
     */
    private String searchForName(int id) {
        try {
            return reverseSearch(this.startOffset + 32, id);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private String reverseSearch(long offset, int id) throws IOException {

        BPTree.Node node = this.bpTree.readTreeNode(offset);

        String found = null;

        if (node.type == 1) {
            // Leaf node
            for (BPTree.Item item : node.items) {
                String key = item.getKey();
                int[] values = item.getValueAsInts();
                nameToId.put(key, values[0]);
                idToName.put(values[0], key);
                if (values[0] == id) {
                    found = key;
                }
            }
            return found;
        } else {
            // non-leaf
            for (BPTree.Item item : node.items) {
                found = reverseSearch(item.offset, id);
                if (found != null) {
                    break;
                }
            }
        }
        return found;
    }


    /**
     * Return an estimated length of the genome, which might be the actual length if the # of contigs is small.
     * This is only used for calculating a default feature visibility window.
     *
     * @return
     */
    public long estimateGenomeSize() {
        try {
            RunningTotal runningTotal = new RunningTotal(0, 0);
            accumulateSize(this.startOffset + 32, runningTotal, 1000);
            return (long) ((double) (this.getItemCount() / runningTotal.count) * runningTotal.total);
        } catch (IOException e) {
            log.error("Error estimating genome size", e);
            return -1;
        }
    }

    private RunningTotal accumulateSize(long offset, RunningTotal runningTotal, int maxCount) throws IOException {

        BPTree.Node node = this.bpTree.readTreeNode(offset);
        if (node.type == 1) {
            // Leaf node
            for (BPTree.Item item : node.items) {
                int[] values = item.getValueAsInts();
                runningTotal.total += values[1];
                runningTotal.count++;
            }
        } else {
            // non-leaf.  Visit items in random order to try to avoid bias.
            List<BPTree.Item> shuffledItems = new ArrayList<>(node.items);
            Collections.shuffle(shuffledItems);
            for (BPTree.Item item : shuffledItems) {
                accumulateSize(item.offset, runningTotal, maxCount);
                if (runningTotal.count > maxCount) {
                    break;
                }
            }
        }
        return runningTotal;
    }
}

class RunningTotal {
    double total;
    int count;

    RunningTotal(double total, int count) {
        this.total = total;
        this.count = count;
    }

    public double getAverage() {
        return this.total / this.count;
    }
}
