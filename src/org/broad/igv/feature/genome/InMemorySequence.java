package org.broad.igv.feature.genome;

import java.io.IOException;
import java.util.Arrays;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: jrobinso
 * Date: 6/23/12
 * Time: 8:59 PM
 * To change this template use File | Settings | File Templates.
 */
public class InMemorySequence {

    private Map<String, byte[]> sequenceMap;

    public InMemorySequence(Map<String, byte[]> sequenceMap) {
        this.sequenceMap = sequenceMap;
    }

    public byte[] readSequence(String chr, int qstart, int qend) {
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

}
