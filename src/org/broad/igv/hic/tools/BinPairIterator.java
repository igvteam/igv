package org.broad.igv.hic.tools;

import org.broad.tribble.util.LittleEndianInputStream;

import java.io.*;
import java.util.Map;
import java.util.zip.GZIPInputStream;

/**
 * @author Jim Robinson
 * @date 4/7/12
 */
public class BinPairIterator implements PairIterator {

    LittleEndianInputStream is;
    AlignmentPair next;

    /**
     * TODO -- chromosomeIndexes is ignored for now, but should be used to map the chromosome stored in the
     * TODO -- bin pair file with an integer index. The current assumption is the chromosome map in
     * TODO -- the bin pair file is the same being used for the hic file,  a fragile assumption.
     *
     * @param path
     * @param chromosomeIndexes
     * @throws IOException
     */
    public BinPairIterator(String path, Map<String, Integer> chromosomeIndexes) throws IOException {
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
        if (is != null) try {
            is.close();
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    private void advance() {

        try {
            boolean str1 = (is.readByte() != 0);
            int chr1 = is.readInt();
            int pos1 = is.readInt();
            int frag1 = is.readInt();
            boolean str2 = (is.readByte() != 0);
            int chr2 = is.readInt();
            int pos2 = is.readInt();
            int frag2 = is.readInt();
            next = new AlignmentPair(str1,chr1, pos1, frag1, str2, chr2, pos2, frag2);
        } catch (IOException e) {
            next = null;
            if (!(e instanceof EOFException)) {
                e.printStackTrace();
            }
        }
    }
}
