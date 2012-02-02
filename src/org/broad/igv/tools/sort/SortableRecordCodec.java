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

package org.broad.igv.tools.sort;

import net.sf.samtools.util.SortingCollection;
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
