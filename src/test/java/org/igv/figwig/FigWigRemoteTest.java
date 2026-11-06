package org.igv.figwig;

import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpServer;
import org.igv.AbstractHeadlessTest;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.net.InetSocketAddress;
import java.nio.file.Files;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicLong;

import static org.junit.Assert.*;

/**
 * figwig files over HTTP range requests, from a server in this JVM: every value and summary equal to
 * the local file's, through redirects and responses shorter than asked; a missing file fails to open.
 */
public class FigWigRemoteTest extends AbstractHeadlessTest {

    static HttpServer server;
    static String base;
    static final AtomicLong requests = new AtomicLong();

    @BeforeClass
    public static void startServer() throws IOException {
        server = HttpServer.create(new InetSocketAddress("127.0.0.1", 0), 64);
        server.createContext("/", FigWigRemoteTest::serve);
        server.setExecutor(Executors.newFixedThreadPool(16));
        server.start();
        base = "http://127.0.0.1:" + server.getAddress().getPort();
    }

    @AfterClass
    public static void stopServer() {
        server.stop(0);
    }

    /** /ok/FILE: ranges; /redirect/FILE: 302 to /ok/FILE; /half/FILE: half of each range asked for. */
    static void serve(HttpExchange x) throws IOException {
        requests.incrementAndGet();
        String[] parts = x.getRequestURI().getPath().split("/", 3);
        File f = parts.length == 3 ? new File(FigWigFileTest.DIR + parts[2]) : null;
        if (f == null || !f.isFile()) {
            x.sendResponseHeaders(404, -1);
            x.close();
            return;
        }
        if (parts[1].equals("redirect")) {
            x.getResponseHeaders().add("Location", "/ok/" + parts[2]);
            x.sendResponseHeaders(302, -1);
            x.close();
            return;
        }
        byte[] all = Files.readAllBytes(f.toPath());
        long start = 0, end = all.length - 1;
        String range = x.getRequestHeaders().getFirst("Range");
        if (range != null) {
            String[] se = range.substring("bytes=".length()).split("-");
            start = Long.parseLong(se[0]);
            end = se.length > 1 && !se[1].isEmpty() ? Math.min(Long.parseLong(se[1]), all.length - 1) : all.length - 1;
            if (start >= all.length) {
                x.getResponseHeaders().add("Content-Range", "bytes */" + all.length);
                x.sendResponseHeaders(416, -1);
                x.close();
                return;
            }
            if (parts[1].equals("half") && end > start) {
                end = start + (end - start + 1) / 2 - 1;
            }
            x.getResponseHeaders().add("Content-Range", "bytes " + start + "-" + end + "/" + all.length);
        }
        x.getResponseHeaders().add("Accept-Ranges", "bytes");
        int n = (int) (end - start + 1);
        if (x.getRequestMethod().equals("HEAD")) {
            x.getResponseHeaders().add("Content-Length", Integer.toString(all.length));
            x.sendResponseHeaders(200, -1);
            x.close();
            return;
        }
        x.sendResponseHeaders(range != null ? 206 : 200, n);
        try (OutputStream out = x.getResponseBody()) {
            out.write(all, (int) start, n);
        }
    }

    static void assertSameAsLocal(String url, String file) throws Exception {
        FigWigFile local = FigWigFile.open(FigWigFileTest.DIR + file);
        FigWigFile remote = FigWigFile.open(url);
        try {
            assertTrue(remote.isHttp());
            assertArrayEquals(local.trackNames, remote.trackNames);
            for (int t = 0; t < local.nTracks(); t++) {
                for (int c = 0; c < local.chromNames.length; c++) {
                    int len = (int) local.chromLengths[c];
                    assertArrayEquals(url + " track " + t + " " + local.chromNames[c],
                            local.readDense(t, c, 0, len), remote.readDense(t, c, 0, len), 0.0);
                    remote.prefetch(t, c, len / 3, len);
                    assertArrayEquals(local.readDense(t, c, len / 3, len), remote.readDense(t, c, len / 3, len), 0.0);
                }
                assertEquals(local.hasSummary(t), remote.hasSummary(t));
                for (int k = 0; k < local.levels(); k++) {
                    FigWigFile.Runs a = local.levelMeans(t, k, 0, 0, (int) local.chromLengths[0]);
                    FigWigFile.Runs b = remote.levelMeans(t, k, 0, 0, (int) local.chromLengths[0]);
                    assertEquals(a == null, b == null);
                    if (a != null) {
                        assertEquals(a.n, b.n);
                        for (int i = 0; i < a.n; i++) {
                            assertEquals(Double.doubleToRawLongBits(a.values[i]), Double.doubleToRawLongBits(b.values[i]));
                        }
                    }
                }
            }
        } finally {
            local.close();
            remote.close();
        }
    }

    @Test(timeout = 300000)
    public void testRangeRequestsReadWhatTheLocalFileHolds() throws Exception {
        for (String file : FigWigFileTest.FILES) {
            assertSameAsLocal(base + "/ok/" + file, file);
        }
    }

    @Test(timeout = 300000)
    public void testRedirectsAndShortResponses() throws Exception {
        for (String file : new String[]{"positions_zlib.fw", "tracks_stored.fw", "cells.fw", "mixed.fw"}) {
            assertSameAsLocal(base + "/redirect/" + file, file);
            assertSameAsLocal(base + "/half/" + file, file);
        }
    }

    @Test(timeout = 60000)
    public void testAMissingFileFailsToOpen() {
        try {
            FigWigFile.open(base + "/ok/missing.fw").close();
            fail("opened a missing file");
        } catch (IOException | RuntimeException expected) {
            // the error names the URL or the HTTP status; nothing hangs
        }
    }
}
