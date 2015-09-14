/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.util;

import org.apache.log4j.Logger;
import org.broad.igv.ui.util.ProgressMonitor;

import java.io.*;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;

/**
 * Runnable for downloading a file from a URL.
 * Downloading is buffered, and can be cancelled (between buffers)
 * via {@link #cancel(boolean)}
 */
public class URLDownloader implements Runnable {

    private static Logger log = Logger.getLogger(URLDownloader.class);

    private ProgressMonitor monitor = null;

    private final URL srcUrl;
    private final File outputFile;

    private volatile boolean started = false;
    private volatile boolean done = false;
    private volatile boolean cancelled = false;
    private volatile RunnableResult result;

    public URLDownloader(String url, File outputFile) throws MalformedURLException {
        this.srcUrl = new URL(url);
        this.outputFile = outputFile;
    }

    @Override
    public void run() {
        if (this.cancelled) {
            return;
        }
        started = true;

        try {
            this.result = doDownload();
        } catch (IOException e) {
            log.error(e.getMessage(), e);
        } finally {
            this.done();
        }

    }

    /**
     * Return the result. Must be called after run is complete
     *
     * @return
     */
    public RunnableResult getResult() {
        if (!this.done) throw new IllegalStateException("Must wait for run to finish before getting result");
        return this.result;
    }

    private RunnableResult doDownload() throws IOException {

       log.info("Downloading " + srcUrl + " to " + outputFile.getAbsolutePath());

        HttpURLConnection conn = HttpUtils.getInstance().openConnection(this.srcUrl, null);

        long contentLength = -1;
        String contentLengthString = conn.getHeaderField("Content-Length");
        if (contentLengthString != null) {
            contentLength = Long.parseLong(contentLengthString);
        }

        InputStream is = null;
        OutputStream out = null;

        long downloaded = 0;
        long downSinceLast = 0;
        String curStatus;
        String msg1 = String.format("downloaded of %s total", contentLength >= 0 ? bytesToByteCountString(contentLength) : "unknown");
        int perc = 0;
        try {
            is = conn.getInputStream();
            out = new FileOutputStream(outputFile);

            byte[] buf = new byte[64 * 1024];
            int counter = 0;
            int interval = 100;
            int bytesRead = 0;
            while (!this.cancelled && (bytesRead = is.read(buf)) != -1) {
                out.write(buf, 0, bytesRead);
                downloaded += bytesRead;
                downSinceLast += bytesRead;
                counter = (counter + 1) % interval;
                if (counter == 0 && this.monitor != null) {
                    curStatus = String.format("%s %s", bytesToByteCountString(downloaded), msg1);
                    this.monitor.updateStatus(curStatus);
                    if (contentLength >= 0) {
                        perc = (int) ((downSinceLast * 100) / contentLength);
                        this.monitor.fireProgressChange(perc);
                        if (perc >= 1) downSinceLast = 0;
                    }
                }
            }
            log.info("Download complete.  Total bytes downloaded = " + downloaded);
        } catch (IOException e) {
            HttpUtils.getInstance().readErrorStream(conn);
            throw e;
        } finally {
            if (is != null) is.close();
            if (out != null) {
                out.flush();
                out.close();
            }
        }
        long fileLength = outputFile.length();

        if (this.cancelled) return RunnableResult.CANCELLED;

        boolean knownComplete = contentLength == fileLength;
        //Assume success if file length not known
        if (knownComplete || contentLength < 0) {
            if (this.monitor != null) {
                this.monitor.fireProgressChange(100);
                this.monitor.updateStatus("Done");
            }
            return RunnableResult.SUCCESS;
        } else {
            return RunnableResult.FAILURE;
        }

    }

    protected void done() {
        this.done = true;
    }

    public boolean isDone() {
        return this.done;
    }

    /**
     * See {@link java.util.concurrent.FutureTask#cancel(boolean)}
     *
     * @param mayInterruptIfRunning
     * @return
     */
    public boolean cancel(boolean mayInterruptIfRunning) {
        if (this.started && !mayInterruptIfRunning) {
            return false;
        }
        this.cancelled = true;
        return true;
    }

    public void setMonitor(ProgressMonitor monitor) {
        this.monitor = monitor;
    }

    /**
     * Convert bytes to human-readable string.
     * e.g. 102894 -> 102.89 KB. If too big or too small,
     * doesn't append a prefix just returns {@code bytes + " B"}
     *
     * @param bytes
     * @return
     */
    public String bytesToByteCountString(long bytes) {
        int unit = 1000;
        String prefs = "KMGT";

        if (bytes < unit) return bytes + " B";
        int exp = (int) (Math.log(bytes) / Math.log(unit));
        if (exp <= 0 || exp >= prefs.length()) return bytes + " B";

        String pre = (prefs).charAt(exp - 1) + "";
        return String.format("%.2f %sB", bytes / Math.pow(unit, exp), pre);
    }
}
