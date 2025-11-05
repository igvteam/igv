package javastraw.reader;

import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.Block;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import org.junit.Test;

import java.util.Iterator;
import java.util.List;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.fail;

public class DatasetReaderV2Test {

    @Test
    public void read() {
        String url = "https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/3be17688-cbce-4ef9-9b94-8571c20a858e/4DNFI916JQ1Y.hic";
        DatasetReaderV2 reader = new DatasetReaderV2(url, false, false);
        try {
            Dataset ds = reader.read();
            assertNotNull(ds);

            List<HiCZoom> zooms = ds.getBpZooms();

            Chromosome c1 = ds.getChromosomeHandler().getChromosomeFromName("1");


            Matrix matrix = ds.getMatrix(c1, c1);

            HiCZoom zoom = zooms.get(5);
            MatrixZoomData zd = matrix.getZoomData(zoom);

            boolean getDataUnderTheDiagonal = false;
            NormalizationType norm = NormalizationHandler.NONE;
            int binXStart = 500, binYStart = 600, binXEnd = 1000, binYEnd = 1200;
            List<Block> blocks = zd.getNormalizedBlocksOverlapping(binXStart, binYStart, binXEnd, binYEnd, norm, getDataUnderTheDiagonal);
            for (Block b : blocks) {
                if (b != null) {
                    for (ContactRecord rec : b.getContactRecords()) {
                        if (rec.getCounts() > 0) { // will skip NaNs
                            // can choose to use the BIN coordinates
                            int binX = rec.getBinX();
                            int binY = rec.getBinY();

                            // you could choose to use relative coordinates for the box given
                            int relativeX = rec.getBinX() - binXStart;
                            int relativeY = rec.getBinY() - binYStart;

                            float counts = rec.getCounts();
                            //System.out.println(binX + "\t" + binY + "\t" + counts);
                        }
                    }
                }
            }

        } catch (Exception e) {
            fail("Exception thrown: " + e.getMessage());
        }
    }

    @Test
    public void readMatrix() {
    }
}