/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam.reader;

import org.broad.igv.util.collections.SortingCollection;
import org.broad.igv.util.collections.SortingCollection.Codec;
import org.broad.igv.logging.*;
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

    private static Logger log = LogManager.getLogger(DotAlignedCodec.class);
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
