package org.igv.sam;

import org.igv.util.TestUtils;
import htsjdk.tribble.readers.AsciiLineReader;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @author jrobinso
 */
public class AsciiLineReaderTest {

    public AsciiLineReaderTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Test
    public void testContent() throws Exception {
        File testFile = new File(TestUtils.DATA_DIR + "igv/recombRate.ens.igv.txt");


        AsciiLineReader reader = new AsciiLineReader(new FileInputStream(testFile));
        BufferedReader br = new BufferedReader(new FileReader(testFile));

        String arLine = null;
        int count = 0;
        while ((arLine = reader.readLine()) != null) {
            String brLine = br.readLine();
            assertEquals(arLine, brLine);
            count++;
        }
        assertTrue(count > 0);
    }

    //AsciiLineReader used to be faster than BufferedReader, no longer.
    //Still may be faster over HTTP
    //TODO Check over HTTP
    @Ignore
    @Test
    public void testSpeed() throws Exception {
        File testFile = new File(TestUtils.DATA_DIR + "cn/HindForGISTIC.hg16.cn");


        AsciiLineReader reader = new AsciiLineReader(new FileInputStream(testFile));

        long asciiCount = 0;
        long t02 = System.currentTimeMillis();
        while (reader.readLine() != null) {
            asciiCount++;
        }
        long asciiReaderTime = System.currentTimeMillis() - t02;

        BufferedReader br = new BufferedReader(new FileReader(testFile));
        long brCount = 0;
        long t0 = System.currentTimeMillis();
        while (br.readLine() != null) {
            brCount++;
        }
        long bufferedReaderTime = System.currentTimeMillis() - t0;


        // It will be considered a bug if AsciiLineReader is slower than BufferedReader
        assertTrue(bufferedReaderTime > asciiReaderTime);
        assertEquals(asciiCount, brCount);

    }
}
