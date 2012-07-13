package org.broad.igv.hic.matrix;

import org.broad.igv.util.ParsingUtils;
import org.broad.tribble.util.LittleEndianInputStream;

import java.io.*;

/**
 * Class for testing/development
 * <p/>
 * Assumptions -- matrix is square
 *
 * @author jrobinso
 *         Date: 7/13/12
 *         Time: 1:48 PM
 */
public class YunfanFormatMatrix implements BasicMatrix {


    String path;
    int dim;
    float[] data;

    // TODO -- store these in file
    String chr = "chr14";
    int binSize = 1000000;


    public YunfanFormatMatrix(float [] data, int dim) {
        this.data = data;
        this.dim = dim;
    }

    public YunfanFormatMatrix(String path) throws IOException {
        this.path = path;
        readMatrix();
    }

    public void readMatrix() throws IOException {
        BufferedInputStream bis = null;
        try {
            InputStream is = ParsingUtils.openInputStream(path);
            bis = new BufferedInputStream(is);
            LittleEndianInputStream les = new LittleEndianInputStream(bis);

            dim = les.readInt();
            int nTot = dim * dim;
            data = new float[nTot];

            for (int i = 0; i < nTot; i++) {
                data[i] = les.readFloat();
                System.out.println(data[i]);
            }

            les.close();
            bis.close();

        } finally {
            if (bis != null)
                bis.close();
        }


    }

    @Override
    public float getEntry(int row, int col) {

        int idx = row * dim + col;
        return data[idx];
    }

    @Override
    public int getRowDimension() {
        return dim;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public int getColumnDimension() {
        return dim;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public BasicMatrix getSubMatrix(int startRow, int endRow, int startCol, int endCol) {

        int startIdx = startRow * dim + startCol;
        int endIdx = endRow * dim + endCol;
        int dim = endRow - startRow + 1;
        // TODO -- assert matrix is square
        float [] subdata = new float[endIdx - startIdx + 1];
        System.arraycopy(data, 0, subdata, 0, subdata.length);

        return new YunfanFormatMatrix(subdata, dim);
    }
}
