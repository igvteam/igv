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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam.reader;

import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.SortingCollection.Codec;
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
