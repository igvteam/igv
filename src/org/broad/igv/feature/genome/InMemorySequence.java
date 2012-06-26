package org.broad.igv.feature.genome;

import java.io.IOException;
import java.util.*;

/**
 * A sequence implementation in which all data is held in memory.
 */
public class InMemorySequence implements Sequence {

    /**
     * Map of chromosome name -> byte array of sequence
     */
    private Map<String, byte[]> sequenceMap;

    public InMemorySequence(Map<String, byte[]> sequenceMap) {
        this.sequenceMap = sequenceMap;
    }

    public InMemorySequence(String chr, byte[] seq) {
        sequenceMap = new HashMap<String, byte[]>();
        sequenceMap.put(chr, seq);
    }

    public byte[] getSequence(String chr, int qstart, int qend) {
        byte[] allBytes = sequenceMap.get(chr);
        if (allBytes == null) {
            return null;
        } else {
            final int start = Math.max(0, qstart);    // qstart should never be < 0
            final int end = Math.min(allBytes.length, qend);
            int len = end - start;

            byte[] bytes = new byte[len];
            Arrays.fill(bytes, (byte) 0);
            int s = Math.max(start, 0);
            System.arraycopy(allBytes, s, bytes, 0, len);
            return bytes;
        }
    }


    @Override
    public byte getBase(String chr, int position) {

        byte[] seqBytes = sequenceMap.get(chr);
        if (seqBytes != null && position < seqBytes.length) {
            return seqBytes[position];
        } else {
            return 0;
        }
    }

    @Override
    public List<String> getChromosomeNames() {
        return new ArrayList<String>(sequenceMap.keySet());
    }

    @Override
    public int getChromosomeLength(String chrname) {
        byte [] bytes = sequenceMap.get(chrname);
        return bytes == null ? 0 : bytes.length;
    }
}
