package org.broad.igv.tools.sort;

import org.broad.igv.logging.*;
import org.broad.igv.util.collections.SortingCollection;

import java.io.*;

/**
 * Codec for Picard sorting classes.   Used to serialize and deserialize records to disk.
 */
public class SortableRecordCodec implements SortingCollection.Codec<SortableRecord> {

    private static Logger log = LogManager.getLogger(SortableRecordCodec.class);
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
