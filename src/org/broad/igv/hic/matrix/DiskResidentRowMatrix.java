package org.broad.igv.hic.matrix;

import org.broad.igv.hic.MainWindow;
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
public class DiskResidentRowMatrix implements BasicMatrix {


    String path;
    int dim;
    float lowerValue;
    float upperValue;
    private String genome;
    private String chr1;
    private String chr2;
    private int binSize;

    int arrayStartPosition;
    boolean isLoading = false;

    ObjectCache<Integer, float[]> rowDataCache = new ObjectCache<Integer, float[]>(10000);

    public DiskResidentRowMatrix(String path) throws IOException {
        this.path = path;
        init();
    }

    public int getDim() {
        return dim;
    }

    public String getGenome() {
        return genome;
    }

    public String getChr1() {
        return chr1;
    }

    public String getChr2() {
        return chr2;
    }

    public int getBinSize() {
        return binSize;
    }

    public float getLowerValue() {
        return lowerValue;
    }

    public float getUpperValue() {
        return upperValue;
    }

    public void loadRowData(int startRow, int lastRow) {
        //System.out.println("Loading row " + row);
        SeekableStream is = null;
        try {
            is = SeekableStreamFactory.getStreamFor(path);

            int startFilePosition = arrayStartPosition + (startRow * dim) * 4;
            int nBytes = (lastRow - startRow) * dim * 4;
            byte[] byteArray = new byte[nBytes];

            is.seek(startFilePosition);
            is.readFully(byteArray);

            ByteArrayInputStream bis = new ByteArrayInputStream(byteArray);
            LittleEndianInputStream les = new LittleEndianInputStream(bis);


            for (int row = startRow; row < lastRow; row++) {
                float[] rowData = new float[dim];
                for (int i = 0; i < dim; i++) {
                    float f = les.readFloat();
                    rowData[i] = f;
                }
                rowDataCache.put(row, rowData);
            }

        } catch (IOException e) {
            e.printStackTrace();
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
    public float getEntry(final int row, int col) {

        if (!rowDataCache.containsKey(row)) {
            MainWindow.getInstance().showGlassPane();
            int lastRow = Math.min(dim, row + 500);
            loadRowData(row, lastRow);
            MainWindow.getInstance().hideGlassPane();
        }
        float[] rowData = rowDataCache.get(row);
        return rowData == null ? Float.NaN : rowData[col];
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

        return new InMemoryMatrix(subdata, dim);
    }


    private void init() throws IOException {
        // Peak at file to determine version
        BufferedInputStream bis = null;
        try {
            InputStream is = ParsingUtils.openInputStream(path);
            bis = new BufferedInputStream(is);
            LittleEndianInputStream les = new LittleEndianInputStream(bis);

            int bytePosition = 0;
            int magic = les.readInt();    // <= 6515048
            bytePosition += 4;
            System.out.println("Magic number = " + magic);

            //if (magic != 6515048)      <= ERROR
            // Version number
            int version = les.readInt();
            bytePosition += 4;
            System.out.println("Version = " + version);

            genome = les.readString();
            bytePosition += genome.length() + 1;
            System.out.println("Genome = " + genome);

            chr1 = les.readString();
            bytePosition += chr1.length() + 1;
            System.out.println("Chr1 = " + chr1);

            chr2 = les.readString();
            bytePosition += chr2.length() + 1;
            System.out.println("Chr2 = " + chr2);

            binSize = les.readInt();
            bytePosition += 4;
            System.out.println("binSize = " + binSize);

            lowerValue = les.readFloat();
            bytePosition += 4;
            System.out.println("Lower value = " + lowerValue);

            upperValue = les.readFloat();
            bytePosition += 4;
            System.out.println("Upper value = " + upperValue);

            int nRows = les.readInt();  // # rows, assuming square matrix
            bytePosition += 4;
            System.out.println("Row count = " + nRows);

            int nCols = les.readInt();
            bytePosition += 4;
            System.out.println("Column count = " + nRows);

            if (nRows == nCols) {
                dim = nRows;
            }
            else {
                throw new RuntimeException("Non-square matrices not supported");
            }

            this.arrayStartPosition = bytePosition;


        } finally {
            if (bis != null)
                bis.close();
        }
    }

}
