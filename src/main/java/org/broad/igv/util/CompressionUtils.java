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

package org.broad.igv.util;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import java.io.*;
import java.util.zip.Deflater;
import java.util.zip.Inflater;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso, mdecautis
 * Date: Dec 13, 2009
 * Time: 4:16:42 PM
 * To change this template use File | Settings | File Templates.
 */
public class CompressionUtils {

    private static Logger log = Logger.getLogger(CompressionUtils.class);

    private Deflater deflater;
    private Inflater decompressor;

    public CompressionUtils() {
        decompressor = new Inflater();
        deflater = new Deflater();
        deflater.setLevel(Deflater.DEFAULT_COMPRESSION);
    }

    public byte[] decompress(byte[] data) {
        return decompress(data, data.length * 4);
    }

    /**
     * @param data                  -- the data to decompress
     * @param uncompressedChunkSize -- an estimate of the uncompressed chunk size.  This need not be exact.
     * @return
     */
    public synchronized byte[] decompress(byte[] data, int uncompressedChunkSize) {

        // mpd: new code
        int rem = data.length;

        // Create an expandable byte array to hold the decompressed data
        ByteArrayOutputStream bos = new ByteArrayOutputStream(uncompressedChunkSize);

        // Decompress the data
        byte[] outbuf = new byte[uncompressedChunkSize];

        decompressor.reset();
        decompressor.setInput(data);
        while (rem > 0) {

            // If we are finished with the current chunk start a new one
            if (decompressor.finished()) {
                decompressor = new Inflater();
                int offset = data.length - rem;
                decompressor.setInput(data, offset, rem);
            }

            try {
                int count = decompressor.inflate(outbuf, 0, outbuf.length);
                rem = decompressor.getRemaining();
                bos.write(outbuf, 0, count);

            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        try {
            bos.close();
        } catch (IOException e) {
            // Ignore -- no resources open
        }

        // Return the decompressed data
        return bos.toByteArray();
    }


    public synchronized byte[] compress(byte[] data) {

        // Give the compressor the data to compress
        deflater.reset();
        deflater.setInput(data);
        deflater.finish();

        // Create an expandable byte array to hold the compressed data.
        // You cannot use an array that's the same size as the orginal because
        // there is no guarantee that the compressed data will be smaller than
        // the uncompressed data.
        ByteArrayOutputStream bos = new ByteArrayOutputStream(data.length);

        // Compress the data
        byte[] buf = new byte[1024];
        while (!deflater.finished()) {
            int count = deflater.deflate(buf);
            bos.write(buf, 0, count);
        }
        try {
            bos.close();
        } catch (IOException e) {
            System.err.println("Error clossing ByteArrayOutputStream");
            e.printStackTrace();
        }

        byte[] compressedData = bos.toByteArray();
        return compressedData;

    }

    public synchronized byte[] compress(byte[] data, int chunkSize) {

        ByteArrayOutputStream bos = new ByteArrayOutputStream(data.length);

        int bytesRemaining = data.length;
        while (bytesRemaining > 0) {
            int sz = Math.min(bytesRemaining, chunkSize);
            int position = data.length - bytesRemaining;
            byte[] chunk = new byte[sz];
            System.arraycopy(data, position, chunk, 0, sz);

            byte[] compressedChunk = compress(chunk);
            bos.write(compressedChunk, 0, compressedChunk.length);

            bytesRemaining -= sz;
        }

        byte[] compressedData = bos.toByteArray();
        return compressedData;

    }

    /**
     * Copy data from srcPath to destPath, ungzipping if necessary.
     * destPath is overwritten if it exists
     *
     * @param srcPath
     * @param destPath If null, strip gz from inputFile
     * @return Where the ungipped file was written to. If destPath != null,
     *         this will be destPath.
     */
    public static String ungzipFile(String srcPath, String destPath) throws IOException {
        InputStream inputStream = null;
        OutputStream outputStream = null;
        if (destPath == null) {
            if (srcPath.endsWith(Globals.GZIP_FILE_EXTENSION)) {
                destPath = srcPath.substring(0, srcPath.length() - Globals.GZIP_FILE_EXTENSION.length());
            } else {
                throw new IllegalArgumentException(srcPath + " does not have a gzip extension and destPath is null. Don't know where to write out");
            }
        }
        try {
            inputStream = ParsingUtils.openInputStream(srcPath);
            outputStream = new FileOutputStream(destPath, false);
            int buffersize = 4096;
            byte[] buffer = new byte[buffersize];

            int len = inputStream.read(buffer);
            while (len > 0) {
                outputStream.write(buffer, 0, len);
                len = inputStream.read(buffer);
            }
            outputStream.flush();
        } catch (IOException e) {
            e.printStackTrace();
            log.error(e.getMessage());
        } finally {
            if (inputStream != null) inputStream.close();
            if (outputStream != null) outputStream.close();
        }
        return destPath;
    }
}
