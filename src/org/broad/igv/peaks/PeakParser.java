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

package org.broad.igv.peaks;


import htsjdk.samtools.seekablestream.SeekableStream;
import org.broad.igv.util.CompressionUtils;
import org.broad.igv.util.stream.IGVSeekableBufferedStream;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;
import htsjdk.tribble.util.LittleEndianInputStream;

import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.IOException;
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
            ss = IGVSeekableStreamFactory.getInstance().getStreamFor(path);
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

                ss = IGVSeekableStreamFactory.getInstance().getStreamFor(path);
                int bufferSize = 512000;
                long contentLength = ss.length();
                if(contentLength > 0) {
                    bufferSize = (int) Math.min(contentLength, bufferSize);
                }

                IGVSeekableBufferedStream bufferedStream = new IGVSeekableBufferedStream(ss, bufferSize);
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

