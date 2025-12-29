package org.igv.util.stream;

import htsjdk.samtools.util.BufferedLineReader;
import org.igv.util.TestUtils;
import org.junit.Test;

import java.io.File;

import static org.junit.Assert.*;

public class IGVSeekableFileStreamTest {


    @Test
    public void testSeek() throws Exception {
        String expectedLine = "ccccccccc";
        File testFile = new File(TestUtils.DATA_DIR + "seekablestream/seekTest.txt");
        IGVSeekableFileStream is = new IGVSeekableFileStream(testFile);
        is.seek(20);
        BufferedLineReader reader = new BufferedLineReader(is);
        String nextLine = reader.readLine();
        assertEquals(expectedLine, nextLine);
        reader.close();
    }
}