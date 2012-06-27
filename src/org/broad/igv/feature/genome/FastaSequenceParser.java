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
public class FastaSequenceParser {


    /**
     * Read an entire fasta file, which might be local or remote and might be gzipped.
     *
     * @param path
     */
    public static Map<String, byte[]> parseFasta(String path) throws IOException {

        Map<String, byte[]> sequenceMap = new HashMap();
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

        return sequenceMap;
    }
}
