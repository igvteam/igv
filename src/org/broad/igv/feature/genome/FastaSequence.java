package org.broad.igv.feature.genome;

import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * A sequence created by reading an entire fasta file (i.e. not an indexed file).
 *
 * @author Jim Robinson
 * @date 3/26/12
 */
public class FastaSequence implements Sequence {

    Map<String, byte[]> sequenceMap;

    public FastaSequence(String path) throws IOException {
        readFasta(path);
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

    /**
     * Read an entire fasta file, which might be local or remote and might be gzipped.
     *
     * @param path
     */
    private void readFasta(String path) throws IOException {

        sequenceMap = new HashMap();
        BufferedReader br = null;

        try {
            br = ParsingUtils.openBufferedReader(path);
            ByteArrayOutputStream buffer = new ByteArrayOutputStream(10000);
            String currentChr = null;
            String nextLine;
            while ((nextLine = br.readLine()) != null) {
                if (nextLine.startsWith("#") || nextLine.trim().length() == 0) {
                    continue;
                } else if (nextLine.startsWith(">")) {
                    if (currentChr != null) {
                        byte[] seq = buffer.toByteArray();
                        sequenceMap.put(currentChr, seq);
                        buffer.reset();   // Resets the count field of this byte array output stream to zero
                    }
                    currentChr = nextLine.substring(1).split("\\s+")[0];
                } else {
                    buffer.write(nextLine.trim().getBytes());
                }
            }
            // Add last chr
            if (currentChr != null) {
                byte[] seq = buffer.toByteArray();
                sequenceMap.put(currentChr, seq);
            }
        } finally {
            if (br != null) br.close();
        }
    }
}
