package org.broad.igv.ucsc.bb;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ucsc.UnsignedByteBuffer;

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

    private Header header;
    private HashMap<String, Integer> nameToId;
    private String[] idToName;

    public long sumLengths = 0;

    public static ChromTree parseTree(UnsignedByteBuffer binaryParser, int startOffset, Genome genome) {

        return (new ChromTree().parse(binaryParser, startOffset, genome));
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

    public ChromTree parse(UnsignedByteBuffer binaryParser, int startOffset, Genome genome) {
        {
            Header header = new Header();
            header.magic = binaryParser.getInt();
            header.blockSize = binaryParser.getInt();
            header.keySize = binaryParser.getInt();
            header.valSize = binaryParser.getInt();
            header.itemCount = binaryParser.getLong();
            header.reserved = binaryParser.getLong();
            this.header = header;

            if (header.valSize != 8) {
                throw new RuntimeException("Unexpected valSize: " + header.valSize);
            }

            this.nameToId = new HashMap<>();
            this.idToName = new String[(int) header.itemCount];

            // Recursively walk tree to populate dictionary
            readTreeNode(-1, startOffset, binaryParser, genome);

            return this;
        }
    }

    void readTreeNode(long offset, int startOffset, UnsignedByteBuffer binaryParser, Genome genome) {

        if (offset >= 0) binaryParser.position((int) offset);
        byte type = binaryParser.get();
        byte reserved = binaryParser.get();
        int count = binaryParser.getUShort();

        if (type == 1) {
            // Leaf node
            for (int i = 0; i < count; i++) {
                String key = binaryParser.getFixedLengthString(header.keySize);
                int value = binaryParser.getInt();
                int chromSize = binaryParser.getInt();
                if (genome != null) {
                    key = genome.getCanonicalChrName(key);  // Translate to canonical chr name
                }
                nameToId.put(key, value);
                idToName[value] = key;
                sumLengths += chromSize;
            }
        } else {
            // non-leaf
            for (int i = 0; i < count; i++) {
                String key = binaryParser.getFixedLengthString(header.keySize);
                long childOffset = binaryParser.getLong();
                long bufferOffset = childOffset - startOffset;
                long currOffset = binaryParser.position();
                readTreeNode(bufferOffset, startOffset, binaryParser, genome);
                binaryParser.position((int) currOffset);
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
