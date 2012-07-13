package org.broad.igv.hic.matrix;

import org.broad.igv.util.ObjectCache;
import org.broad.igv.util.ParsingUtils;
import org.broad.tribble.util.LittleEndianInputStream;
import org.broad.tribble.util.SeekableStream;
import org.broad.tribble.util.SeekableStreamFactory;

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
public class YunfanFormatMatrix2 implements BasicMatrix {


    String path;
    int dim;

    ObjectCache<Integer, float[]> rowDataCache = new ObjectCache<Integer, float[]>(10000);

    // TODO -- store these in file
    String chr = "chr14";
    int binSize = 1000000;


    public YunfanFormatMatrix2(String path) throws IOException {
        this.path = path;
        init();
    }

    public void init() throws IOException {
        BufferedInputStream bis = null;
        try {
            InputStream is = ParsingUtils.openInputStream(path);
            bis = new BufferedInputStream(is);
            LittleEndianInputStream les = new LittleEndianInputStream(bis);

            dim = les.readInt();

            les.close();
            bis.close();

        } finally {
            if (bis != null)
                bis.close();
        }
    }

    public float[] loadRowData(int row) {
        System.out.println("Loading row " + row);
        SeekableStream is = null;
        try {
            is = SeekableStreamFactory.getStreamFor(path);

            int startFilePosition = 4 * 1 + (row * dim) * 4;
            is.seek(startFilePosition);
            BufferedInputStream bis = new BufferedInputStream(is);
            LittleEndianInputStream les = new LittleEndianInputStream(bis);

            float[] rowData = new float[dim];

            for (int i = 0; i < dim; i++) {
                rowData[i] = les.readFloat();
            }

            les.close();
            bis.close();

            return rowData;

        } catch (IOException e) {
            e.printStackTrace();
            return null;
        } finally {
            if (is != null)
                try {
                    is.close();
                } catch (IOException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
        }
    }

    @Override
    public float getEntry(int row, int col) {

        float[] rowData = rowDataCache.get(row);
        if (rowData == null) {
            rowData = loadRowData(row);
            rowDataCache.put(row, rowData);
        }

        return rowData[col];
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
        float[] subdata = new float[endIdx - startIdx + 1];

        for (int row = startRow; row <= endRow; row++) {
            for (int col = startCol; col <= endCol; col++) {
                int idx = (row - startRow) * dim + row + (col - startCol) + col;
                subdata[idx] = getEntry(row, col);
            }
        }

        return new YunfanFormatMatrix(subdata, dim);
    }
}
