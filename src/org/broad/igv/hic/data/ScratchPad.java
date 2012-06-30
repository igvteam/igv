package org.broad.igv.hic.data;

import org.broad.tribble.util.LittleEndianInputStream;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

/**
 * @author jrobinso
 *         Date: 6/29/12
 *         Time: 2:27 PM
 */
public class ScratchPad {

    public static void main(String [] args) throws IOException {
        File f = new File("/Users/jrobinso/bin_chr14_1M.bin");
        readData(f);
    }

    public static void readData(File file ) throws IOException {

        FileInputStream fis = null;

        fis = new FileInputStream(file);
        BufferedInputStream bis = new BufferedInputStream(fis);
        LittleEndianInputStream les = new LittleEndianInputStream(bis);

        int nRows = les.readInt();
        int nTot = nRows*nRows;
        int nBytes = 4;
        for(int i=0; i<nTot; i++) {
            //les.readByte();
            Float f = les.readFloat();
            System.out.println(f);
            nBytes += 4;

        }

        System.out.println("Nbytes=" + nBytes);


        fis.close();
    }
}
