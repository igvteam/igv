package org.igv.util.stream;

import junit.framework.TestCase;
import org.igv.util.HttpUtils;
import org.junit.Test;

import java.net.URL;

/**
 * @author jrobinso
 * @date Jul 28, 2010
 */
public class SeekableServiceStreamTest extends TestCase {

    /**
     * Test a file at some random position using the webservice, and compare results obtained to the standard
     * http stream.
     *
     * @throws Exception
     */
    @Test
    public void testRead() throws Exception {

        String testFile = "https://1000genomes.s3.amazonaws.com/phase3/data/HG01894/exome_alignment/HG01894.mapped.ILLUMINA.bwa.ACB.exome.20130415.bam";

        //HttpUtils.getInstance().updateProxySettings();

        IGVSeekableHTTPStream hs = new IGVSeekableHTTPStream(HttpUtils.createURL(testFile));
        final int position = 100;
        hs.seek(position);
        final int range = 1000;
        byte[] expectedBytes = new byte[range];
        hs.read(expectedBytes, 0, expectedBytes.length);

        SeekableServiceStream sss = new SeekableServiceStream(HttpUtils.createURL(testFile));
        sss.seek(position);
        byte[] bytes = new byte[range];
        sss.read(bytes, 0, bytes.length);

        for (int i = 0; i < expectedBytes.length; i++) {
            assertEquals(expectedBytes[i], bytes[i]);
        }
    }

    /**
     * Test reading past the end of the file.  The service should return bytes up to end of the file.
     *
     * @throws Exception
     */
    @Test
    public void testReadOver() throws Exception {

        String testFile = "https://1000genomes.s3.amazonaws.com/phase3/data/HG01894/exome_alignment/HG01894.mapped.ILLUMINA.bwa.ACB.exome.20130415.bam";

        HttpUtils.getInstance().updateProxySettings();

        long contentLength = HttpUtils.getInstance().getContentLength(new URL(testFile));

        IGVSeekableHTTPStream hs = new IGVSeekableHTTPStream(HttpUtils.createURL(testFile));
        final long position = contentLength - 500;
        hs.seek(position);
        final int range = 1000;
        byte[] expectedBytes = new byte[range];
        int bytesRead = hs.read(expectedBytes, 0, expectedBytes.length);

        SeekableServiceStream sss = new SeekableServiceStream(HttpUtils.createURL(testFile));
        sss.seek(position);
        byte[] bytes = new byte[range];
        int bytesRead2 = sss.read(bytes, 0, bytes.length);

        assertEquals(500, bytesRead);
        assertEquals(500, bytesRead2);
        for (int i = 0; i < expectedBytes.length; i++) {
            assertEquals(expectedBytes[i], bytes[i]);
        }
    }


}
