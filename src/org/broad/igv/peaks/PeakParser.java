/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.peaks;


import org.broad.igv.util.CompressionUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.util.LittleEndianInputStream;
import org.broad.tribble.util.SeekableBufferedStream;
import org.broad.tribble.util.SeekableStream;
import org.broad.tribble.util.SeekableStreamFactory;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 * @date Apr 22, 2011
 */
public class PeakParser {

    Map<String, Long> chrIndex;
    int nTimePoints;
    String path;
    String trackLine;
    int[] times;
    String signalsPath;
    String[] timeSignalsPath;
    private final CompressionUtils compressionUtils;

    public PeakParser(String path) throws IOException {
        this.path = path;
        loadHeader();
        compressionUtils = new CompressionUtils();
    }

    private void loadHeader() throws IOException {
        SeekableStream ss = null;
        try {
            ss = SeekableStreamFactory.getStreamFor(path);
            LittleEndianInputStream is = new LittleEndianInputStream(new BufferedInputStream(ss));

            long indexPosition = is.readLong();
            trackLine = is.readString();
            nTimePoints = is.readInt();
            times = new int[nTimePoints];
            for (int t = 0; t < nTimePoints; t++) {
                times[t] = is.readInt();
            }
            signalsPath = is.readString();
            timeSignalsPath = new String[nTimePoints];
            for (int t = 0; t < nTimePoints; t++) {
                timeSignalsPath[t] = is.readString();
            }

            chrIndex = new HashMap<String, Long>();

            ss.seek(indexPosition);
            is = new LittleEndianInputStream(new BufferedInputStream(ss));
            int nChrs = is.readInt();
            for (int i = 0; i < nChrs; i++) {
                String chr = is.readString();
                long pos = is.readLong();
                chrIndex.put(chr, pos);
            }

        } finally {
            if (ss != null) ss.close();
        }

    }


    public List<Peak> loadPeaks(String chr) throws IOException {


        Long chrPos = chrIndex.get(chr);
        if (chrPos == null) {
            return new ArrayList<Peak>();
        } else {
            List<Peak> peaks = new ArrayList<Peak>(10000);
            LittleEndianInputStream reader = null;
            SeekableStream ss = null;
            try {

                ss = SeekableStreamFactory.getStreamFor(path);
                int bufferSize = 512000;
                long contentLength = ss.length();
                if(contentLength > 0) {
                    bufferSize = (int) Math.min(contentLength, bufferSize);
                }

                SeekableBufferedStream bufferedStream = new SeekableBufferedStream(ss, bufferSize);
                bufferedStream.seek(chrPos);

                reader = new LittleEndianInputStream(bufferedStream);
                int nBytes = reader.readInt();

                byte[] compressedBytes = new byte[nBytes];
                bufferedStream.readFully(compressedBytes);

                byte[] bytes = compressionUtils.decompress(compressedBytes);

                ByteArrayInputStream bis = new ByteArrayInputStream(bytes);
                reader = new LittleEndianInputStream(bis);

                String chrRecorded = reader.readString();
                if (!chrRecorded.equals(chr)) {
                    throw new RuntimeException("Error paring peak file: " + path +
                            "<br>Expected: " + chr + "  found: " + chrRecorded);
                }

                int nDataPoints = reader.readInt();
                for (int n = 0; n < nDataPoints; n++) {
                    int start = reader.readInt();
                    int end = reader.readInt();
                    float combinedScore = reader.readFloat();
                    float[] timePointScores = new float[nTimePoints];
                    for (int i = 0; i < nTimePoints; i++) {
                        timePointScores[i] = reader.readFloat();
                    }
                    peaks.add(new Peak(chr, start, end, "", combinedScore, timePointScores));

                }
                return peaks;

            } finally {
                if (ss != null) ss.close();
            }
        }

    }

}

