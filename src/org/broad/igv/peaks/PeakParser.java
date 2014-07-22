/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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

