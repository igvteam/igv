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
