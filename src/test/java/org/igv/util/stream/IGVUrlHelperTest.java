package org.igv.util.stream;

import com.sun.net.httpserver.HttpServer;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureReader;
import org.igv.AbstractHeadlessTest;
import org.igv.track.TribbleFeatureSource;
import org.igv.util.ResourceLocator;
import org.junit.Test;

import java.io.ByteArrayOutputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.net.InetSocketAddress;

import static org.junit.Assert.assertEquals;

public class IGVUrlHelperTest extends AbstractHeadlessTest {

    /**
     * Gzipped (non-indexed) feature files are decoded by htsjdk with a GZIPInputStream, which before JDK 23 can stop
     * at a bgzf block boundary when the network stream momentarily has no bytes available.
     * See https://github.com/igvteam/igv/issues/1693
     */
    @Test
    public void testRemoteBgzippedBed() throws Exception {

        int nFeatures = 5000;   // Several 64kb bgzf blocks
        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        try (PrintWriter pw = new PrintWriter(new BlockCompressedOutputStream(bos, (java.nio.file.Path) null))) {
            for (int i = 0; i < nFeatures; i++) {
                pw.println("chr1\t" + (i * 100) + "\t" + (i * 100 + 50) + "\tfeature" + i);
            }
        }
        byte[] bytes = bos.toByteArray();

        // Serve the file one bgzf block at a time, pausing between blocks so available() is 0 at block boundaries.
        HttpServer server = HttpServer.create(new InetSocketAddress("localhost", 0), 0);
        server.createContext("/", exchange -> {
            if (!exchange.getRequestURI().getPath().endsWith(".bed.gz")) {
                exchange.sendResponseHeaders(404, -1);
            } else if (exchange.getRequestMethod().equals("HEAD")) {
                exchange.sendResponseHeaders(200, -1);
            } else {
                exchange.sendResponseHeaders(200, 0);
                try (OutputStream out = exchange.getResponseBody()) {
                    int offset = 0;
                    while (offset < bytes.length) {
                        int blockSize = ((bytes[offset + 16] & 0xff) | ((bytes[offset + 17] & 0xff) << 8)) + 1;
                        out.write(bytes, offset, blockSize);
                        out.flush();
                        offset += blockSize;
                        Thread.sleep(300);
                    }
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                }
            }
            exchange.close();
        });
        server.start();

        htsjdk.tribble.util.ParsingUtils.setURLHelperFactory(IGVUrlHelperFactory.getInstance());
        try {
            String url = "http://localhost:" + server.getAddress().getPort() + "/multiblock.bed.gz";
            FeatureReader<Feature> reader =TribbleFeatureSource.getBasicReader(new ResourceLocator(url), genome);
            int count = 0;
            try (CloseableTribbleIterator<Feature> iter = reader.iterator()) {
                while (iter.hasNext()) {
                    iter.next();
                    count++;
                }
            }
            assertEquals(nFeatures, count);
        } finally {
            server.stop(0);
        }
    }
}
