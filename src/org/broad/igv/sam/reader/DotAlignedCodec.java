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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam.reader;

import net.sf.samtools.util.SortingCollection;
import net.sf.samtools.util.SortingCollection.Codec;
import org.apache.log4j.Logger;
import org.broad.igv.sam.DotAlignedAlignment;

import java.io.*;

/**
 * private String chromosome;
 * private int start;
 * private int end;
 * boolean negativeStrand;
 *
 * @author jrobinso
 */
public class DotAlignedCodec implements SortingCollection.Codec<DotAlignedAlignment> {

    private static Logger log = Logger.getLogger(DotAlignedCodec.class);
    DataOutputStream outputStream;
    DataInputStream inputStream;

    public void setOutputStream(OutputStream outputStream) {
        this.outputStream = new DataOutputStream(outputStream);
    }

    public void setInputStream(InputStream inputStream) {
        this.inputStream = new DataInputStream(inputStream);
    }

    public void encode(DotAlignedAlignment alignment) {
        try {
            outputStream.writeUTF(alignment.getChr());
            outputStream.writeInt(alignment.getStart());
            outputStream.writeInt(alignment.getEnd());
            outputStream.writeBoolean(alignment.isNegativeStrand());
        } catch (IOException ex) {
            log.error("Error encoding alignment", ex);
        }
    }

    public DotAlignedAlignment decode() {
        try {

            String chr = inputStream.readUTF();
            int start = inputStream.readInt();
            int end = inputStream.readInt();
            boolean negativeStrand = inputStream.readBoolean();
            return new DotAlignedAlignment(chr, start, end, negativeStrand);
        } catch (EOFException ex) {
            return null;
        } catch (IOException ex) {
            log.error("Error decoding alignment", ex);
            return null;
        }
    }

    public Codec<DotAlignedAlignment> clone() {
        DotAlignedCodec other = new DotAlignedCodec();
        return other;
    }
}
