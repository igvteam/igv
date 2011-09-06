package org.broad.igv.hic.data;

import net.sf.samtools.util.SortingCollection;

import java.io.*;

/**
 * Codec provide for sorting utility, not necessarily for permanent storage
 *
 * @author jrobinso
 * @date Aug 4, 2010
 */
public class ContactRecordCodec implements SortingCollection.Codec<ContactRecord> {

    DataOutputStream os;
    DataInputStream is;



    public void setOutputStream(OutputStream os) {
       this.os = new DataOutputStream(os);
    }

    public void setInputStream(InputStream is) {
        this.is = new DataInputStream(is);
    }

    public void encode(ContactRecord contactRecord) {
        try {
            os.writeInt(contactRecord.getBlockNumber());
            os.writeInt(contactRecord.getX());
            os.writeInt(contactRecord.getY());
            os.writeShort(contactRecord.getCounts());
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

    }

    public ContactRecord decode() {
        try {
            int block = is.readInt();
            int bin1 = is.readInt();
            int bin2 = is.readInt();
            short counts = is.readShort();
            return new ContactRecord(block, bin1, bin2, counts);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            throw new RuntimeException(e);
        }

    }

    public SortingCollection.Codec<ContactRecord> clone() {
        ContactRecordCodec newCodec = new ContactRecordCodec();
        newCodec.os = os;
        newCodec.is = is;
        return newCodec;
    }
}
