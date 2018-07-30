package org.broad.igv.ui.util.download;


import org.apache.log4j.Logger;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.HttpUtils;

import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.net.*;
import java.nio.file.CopyOption;
import java.nio.file.Files;
import java.util.*;

import static java.nio.file.StandardCopyOption.REPLACE_EXISTING;


// This class downloads a file from a URL.
public class Downloader implements Runnable {

    private static final Logger log = Logger.getLogger(Downloader.class);

    // Max size of download buffer.
    private static final int MAX_BUFFER_SIZE = 1000000;    // Max buffer size


    private URL url; // download URL
    private File localFile;
    private String tmpName;
    private int downloaded; // number of bytes downloaded
    private boolean canceled;
    private final ProgressMonitor monitor;

    // Constructor for Download.
    public Downloader(URL url, File file, ProgressMonitor monitor) {
        this.url = url;
        this.localFile = file;
        this.tmpName = file + ".download";
        this.monitor = monitor;
        this.canceled = false;
        this.downloaded = 0;
    }


    // Download file.
    public void run() {

        RandomAccessFile file = null;
        InputStream stream = null;

        try {

            long contentLength = HttpUtils.getInstance().getContentLength(url);


            // Check for valid content length.
            if (contentLength < 1) {
                // error();
            }

            Map<String, String> requestProperties = new HashMap<>();
            requestProperties.put("Range", "bytes=" + downloaded + "-" + downloaded + MAX_BUFFER_SIZE);

            // Open connection to URL.
            HttpURLConnection connection = HttpUtils.getInstance().openConnection(url, null);


            // Connect to server.
            connection.connect();

            // Make sure response code is in the 200 range.
            final int responseCode = connection.getResponseCode();
            if (responseCode < 200 || responseCode >= 300) {
                MessageUtils.showMessage("Error downloading " + url + ". Response code: " + responseCode);
                return;
            }


            // Open file and seek to the end of it.
            file = new RandomAccessFile(tmpName, "rw");
            file.seek(downloaded);

            stream = connection.getInputStream();
            while (downloaded < contentLength) {

                // Size buffer according to how much of the file is left to download.
                byte buffer[];
                if (contentLength - downloaded > MAX_BUFFER_SIZE) {
                    buffer = new byte[MAX_BUFFER_SIZE];
                } else {
                    buffer = new byte[(int) (contentLength - downloaded)];
                }

                // Read from server into buffer.
                int read = stream.read(buffer);
                if (read == -1)
                    break;

                // Write buffer to file.
                file.write(buffer, 0, read);
                downloaded += read;

                int percent = (int) (100.0 * downloaded / contentLength);

                if (monitor != null) {

                    if (monitor.isCanceled()) {
                        this.canceled = true;
                        break;
                    } else {
                        SwingUtilities.invokeLater(() -> {
                            monitor.setProgress(percent);
                            monitor.setNote("" + (downloaded / 1000) + " of " + (contentLength / 1000) + " kb");
                        });
                    }
                }
            }
            file.close();

        } catch (Exception e) {
            log.error("Error downloading " + url, e);

        } finally {

            if (monitor != null) {
                SwingUtilities.invokeLater(() -> monitor.close());
   //             monitor.close();
            }

            // Close file.
            if (file != null) {
                try {
                    file.close();
                } catch (Exception e) {
                }
            }

            // Close connection to server.
            if (stream != null) {
                try {
                    stream.close();
                } catch (Exception e) {
                }
            }


            if(canceled) {
                (new File(tmpName)).delete();
            }
            else {
                try {
                    Files.move(new File(tmpName).toPath(), localFile.toPath(), REPLACE_EXISTING);
                } catch (IOException e) {
                    log.error("Error renaming download file " + tmpName, e);
                }
            }

        }

    }

    /**
     * Convenience method.  Download the resource at url to the local file.
     *
     * @param url
     * @param localFile
     * @param frame
     * @return  true for completed download, false if canceled.
     * @throws MalformedURLException
     */
    public static boolean download(URL url, File localFile, Component frame) throws MalformedURLException {

        String message = "Downloading " + url.toString();
        int min = 0;
        int max = 100;

        final javax.swing.ProgressMonitor monitor;
        if(IGV.hasInstance()) {
            monitor = new javax.swing.ProgressMonitor(frame, message, "", min, max);
            monitor.setMillisToDecideToPopup(100);
        } else {
            monitor = null;
        }

        final Downloader downloader = new Downloader(url, localFile, monitor);
        downloader.run();

        return downloader.canceled == false;


    }

    public static void main(String[] args) throws MalformedURLException {

        URL url = HttpUtils.createURL(args[0]);
        String localFile = args[1];
        JComponent frame = null;

        String message = "Downloading " + url.toString();
        int min = 0;
        int max = 100;

        final javax.swing.ProgressMonitor monitor = new javax.swing.ProgressMonitor(frame, message, "", min, max);
        monitor.setMillisToDecideToPopup(100);


        Downloader dl = new Downloader(url, new File(localFile), monitor);

        (new Thread(dl)).start();


    }

}
