package org.broad.igv.hic.tools;

import org.broad.tribble.util.LittleEndianInputStream;

import java.io.*;
import java.util.zip.GZIPInputStream;

/**
 * @author Jim Robinson
 * @date 4/7/12
 */
public class BinPairIterator implements PairIterator {

    LittleEndianInputStream is;
    AlignmentPair next;

    public BinPairIterator(String path) throws IOException {
        is = new LittleEndianInputStream(new BufferedInputStream(new FileInputStream(path)));
        advance();
    }

    public boolean hasNext() {
        return next != null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public AlignmentPair next() {
        AlignmentPair retValue = next;
        advance();
        return retValue;
    }

    public void remove() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void close() {
        if(is != null) try {
            is.close();
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    private void advance() {

        try {
            int chr1 = is.readInt();
            int pos1 = is.readInt();
            int chr2 = is.readInt();
            int pos2 = is.readInt();
            next = new AlignmentPair(chr1, pos1, chr2, pos2);
        } catch (IOException e) {
            next = null;
            if(!(e instanceof EOFException)) {
                e.printStackTrace();
            }
        }
    }
}
