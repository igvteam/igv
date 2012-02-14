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
package org.broad.igv.bbfile;
import org.apache.log4j.Logger;

import java.util.zip.Inflater;
import java.util.zip.Deflater;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: martind
 * Date: Mar 11, 2010
 * Time: 11:00:57 AM
 * To change this template use File | Settings | File Templates.
 */
/*
*   Compression/Decompression Utillity adapted to BigBed/BigWig Compression formats.
*
*   Note: modified from igv CompressionUtil
*
* */
public class BBCompressionUtils {

    private static Logger log = Logger.getLogger(BBCompressionUtils.class);

    /*
    * Decompress ZLIB commpressed data into a buffer
    *
    * Parameters:
    *   data - data buffer, where data.length indicates requested number of
    *       decompressed bytes.
    *   uncompressBufSize - recommended size for decompression chunks, or can
    *       be full byte size for the compressed region
    *
    * Return:
    *   buffer of uncompressed byte data
    * */
    public static byte[] decompress(byte[] data, int uncompressBufSize) {

        // mpd: new code
        byte [] inbuf = data;   // input is first assigned to full data
        int rem = data.length;
        int count = 0;
        int off = 0;

        // Create an expandable byte array to hold the decompressed data
        ByteArrayOutputStream bos = new ByteArrayOutputStream((int) (1.5 * inbuf.length));

        // Decompress the data
        byte[] outbuf = new byte[uncompressBufSize];

        while (rem > 0) {

            // Create the decompressor and give it the data to compress
            Inflater decompressor = new Inflater();
            decompressor.setInput(inbuf);

            try {
                //byte[] outbuf = new byte[uncompressBufSize];
                uncompressBufSize = outbuf.length;
                count = decompressor.inflate(outbuf, 0, uncompressBufSize);

                rem = decompressor.getRemaining();
                if (rem > 0) {
                    off += inbuf.length - rem;
                    inbuf = new byte[rem];
                    System.arraycopy(data, off, inbuf, 0, rem);
                    uncompressBufSize -= count;
                }

                bos.write(outbuf, 0, count);

            //} catch (DataFormatException e) {
            } catch (Exception e) {
                log.error(e.getMessage());
            }
        }

        try {
            bos.close();
        } catch (IOException e) {
        }

        // Get the decompressed data
        return bos.toByteArray();
    }
    /*
    * */
    public static byte[] compress(byte[] data, int compressBufSize) {
        Deflater compressor = new Deflater();
        compressor.setLevel(Deflater.DEFAULT_COMPRESSION);

        // Give the compressor the data to compress
        compressor.setInput(data);
        compressor.finish();

        // Create an expandable byte array to hold the compressed data.
        // You cannot use an array that's the same size as the orginal because
        // there is no guarantee that the compressed data will be smaller than
        // the uncompressed data.
        ByteArrayOutputStream bos = new ByteArrayOutputStream(data.length);

        // Compress the data
        byte[] buf = new byte[compressBufSize];
        while (!compressor.finished()) {
            int count = compressor.deflate(buf);
            bos.write(buf, 0, count);
        }
        try {
            bos.close();
        } catch (IOException e) {
        }

        byte[] compressedData = bos.toByteArray();
        return compressedData;

    }

    
}
