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

package org.broad.igv.tools.sort;

import htsjdk.samtools.util.SortingCollection;
import org.apache.log4j.Logger;

import java.io.*;

/**
 * Codec for Picard sorting classes.   Used to serialize and deserialize records to disk.
 */
public class SortableRecordCodec implements SortingCollection.Codec<SortableRecord> {

    private static Logger log = Logger.getLogger(SortableRecordCodec.class);
    DataOutputStream outputStream;
    DataInputStream inputStream;

    public void setOutputStream(OutputStream outputStream) {
        this.outputStream = new DataOutputStream(outputStream);
    }

    public void setInputStream(InputStream inputStream) {
        this.inputStream = new DataInputStream(inputStream);
    }

    public void encode(SortableRecord record) {
        try {
            outputStream.writeUTF(record.getChromosome());
            outputStream.writeInt(record.getStart());

            // Code below contributed by Eric Smith to deal with lines > 64k characters (VCF files), which causes
            // writeUTF to blow up since it uses 16-bit length.  The workaround is to write a 32-bit
            // length followed by the UTF8-encoded bytes.
            String s = record.getText();
            byte[] textBytes = s.getBytes("utf-8");
            outputStream.writeInt(textBytes.length);
            outputStream.write(textBytes, 0, textBytes.length);
        } catch (IOException ex) {
            log.error("Error encoding alignment", ex);
        }
    }

    public SortableRecord decode() {
        try {

            String chr = inputStream.readUTF();
            int start = inputStream.readInt();

            // See comment in encode re long lines and writeUTF
            int textLen = inputStream.readInt();
            byte[] textBytes = new byte[textLen];
            inputStream.readFully(textBytes);
            String text = new String(textBytes, "utf-8");

            return new SortableRecord(chr, start, text);
        } catch (EOFException ex) {
            return null;
        } catch (IOException ex) {
            log.error("Error decoding alignment", ex);
            return null;
        }
    }

    public SortingCollection.Codec<SortableRecord> clone() {
        SortableRecordCodec other = new SortableRecordCodec();
        return other;
    }
}
