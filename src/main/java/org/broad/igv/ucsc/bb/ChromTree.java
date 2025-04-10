package org.broad.igv.ucsc.bb;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ucsc.twobit.UnsignedByteBuffer;

import java.util.HashMap;

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

    Header header;
    private HashMap<String, Integer> nameToId;
    private String[] idToName;

    public long sumLengths = 0;

    public static ChromTree parseTree(UnsignedByteBuffer buffer, long startOffset, Genome genome) {

        return (new ChromTree().parse(buffer, startOffset, genome));
    }

    private ChromTree() {

    }

    public Integer getIdForName(String chr) {
        return nameToId.get(chr);
    }

    public String getNameForId(int id) {
        if (id < 0 || id >= idToName.length) {
            return null;
        } else {
            return idToName[id];
        }
    }

    public String[] names() {
        return idToName;
    }

    public ChromTree parse(UnsignedByteBuffer buffer, long startOffset, Genome genome) {
        {
            Header header = new Header();
            header.magic = buffer.getInt();
            header.blockSize = buffer.getInt();
            header.keySize = buffer.getInt();
            header.valSize = buffer.getInt();
            header.itemCount = buffer.getLong();
            header.reserved = buffer.getLong();
            this.header = header;

            if (header.valSize != 8) {
                throw new RuntimeException("Unexpected valSize: " + header.valSize);
            }

            this.nameToId = new HashMap<>();
            this.idToName = new String[(int) header.itemCount];

            // Recursively walk tree to populate dictionary
            readTreeNode(-1, startOffset, buffer, genome);

            return this;
        }
    }

    void readTreeNode(long offset, final long startOffset, UnsignedByteBuffer buffer, Genome genome) {
        if (offset >= 0) buffer.position((int) offset);
        byte type = buffer.get();
        byte reserved = buffer.get();
        int count = buffer.getUShort();

        if (type == 1) {
            // Leaf node
            for (int i = 0; i < count; i++) {
                String key = buffer.getFixedLengthString(header.keySize);
                int value = buffer.getInt();
                int chromSize = buffer.getInt();
                nameToId.put(key, value);
                idToName[value] = key;
                sumLengths += chromSize;
            }
        } else {
            // non-leaf
            for (int i = 0; i < count; i++) {
                String key = buffer.getFixedLengthString(header.keySize);
                long childOffset = buffer.getLong();
                long bufferOffset = childOffset - startOffset;
                long currOffset = buffer.position();
                readTreeNode(bufferOffset, startOffset, buffer, genome);
                buffer.position((int) currOffset);
            }
        }
    }


    static class Header {
        int magic;
        int blockSize;
        int keySize;
        int valSize;
        long itemCount;
        long reserved;
    }
}
