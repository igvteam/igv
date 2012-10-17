/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.sam.reader;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.SeekableBufferedStream;
import net.sf.samtools.util.SeekableStream;
import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.sam.Alignment;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;
import org.broad.igv.util.stream.SeekablePicardStream;
import org.broad.tribble.util.SeekableFTPStream;

import java.io.*;
import java.net.URL;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Sep 22, 2009
 * Time: 2:21:04 PM
 * To change this template use File | Settings | File Templates.
 */
public class BAMHttpReader implements AlignmentReader {

    static Logger log = Logger.getLogger(BAMHttpReader.class);

    // Length of day in milliseconds
    public static final long oneDay = 24 * 60 * 60 * 1000;

    static Hashtable<String, File> indexFileCache = new Hashtable<String, File>();

    URL url;
    SAMFileHeader header;
    File indexFile;
    SAMFileReader reader;
    List<String> sequenceNames;

    public BAMHttpReader(ResourceLocator locator, boolean requireIndex) throws IOException {
        this.url = new URL(locator.getPath());
        if (requireIndex) {
            indexFile = getIndexFile(url, locator.getIndexPath());
            if (indexFile == null) {
                throw new RuntimeException("Could not load index file for file: " + url.getPath());
            }
            //SeekableStream ss = new SeekableBufferedStream(getSeekableStream(url));
            SeekableStream ss = getSeekableStream(url);
            reader = new SAMFileReader(ss, indexFile, false);
        } else {
            InputStream is = HttpUtils.getInstance().openConnectionStream(url);
            reader = new SAMFileReader(new BufferedInputStream(is));
        }

    }

    public void close() throws IOException {
        if (reader != null) {
            reader.close();
        }
    }

    public SAMFileHeader getHeader() {
        if (header == null) {
            header = reader.getFileHeader();
        }
        return header;
    }

    public Set<String> getPlatforms() {
        return AlignmentReaderFactory.getPlatforms(getHeader());
    }

    public boolean hasIndex() {
        return indexFile != null && indexFile.exists();
    }

    public List<String> getSequenceNames() {
        if (sequenceNames == null) {
            SAMFileHeader header = getHeader();
            if (header == null) {
                return null;
            }
            sequenceNames = new ArrayList();
            List<SAMSequenceRecord> records = header.getSequenceDictionary().getSequences();
            if (records.size() > 0) {
                for (SAMSequenceRecord rec : header.getSequenceDictionary().getSequences()) {
                    String chr = rec.getSequenceName();
                    sequenceNames.add(chr);
                }
            }
        }
        return sequenceNames;
    }


    public CloseableIterator<Alignment> iterator() {
        try {
            if (reader == null) {
                InputStream is = HttpUtils.getInstance().openConnectionStream(url);
                reader = new SAMFileReader(new BufferedInputStream(is));
            }
            return new WrappedIterator(reader.iterator());
        } catch (IOException e) {
            log.error("Error creating iterator", e);
            throw new RuntimeException(e);
        }

    }

    public CloseableIterator<Alignment> query(String sequence, int start, int end, boolean contained) {
        try {
            if (reader == null) {
                SeekableStream ss = new SeekableBufferedStream(getSeekableStream(url));
                reader = new SAMFileReader(ss, indexFile, false);
            }
            CloseableIterator<SAMRecord> iter = reader.query(sequence, start + 1, end, contained);
            return new WrappedIterator(iter);
        } catch (IOException e) {
            log.error("Error opening SAM reader", e);
            throw new RuntimeException("Error opening SAM reader", e);
        }
    }

    private SeekableStream getSeekableStream(URL url) throws IOException {
        String protocol = url.getProtocol().toLowerCase();
        SeekableStream is = null;
        if (protocol.equals("http") || protocol.equals("https")) {
            org.broad.tribble.util.SeekableStream tribbleStream = IGVSeekableStreamFactory.getStreamFor(url.toExternalForm());
            String source = url.toExternalForm();
            is = new SeekablePicardStream(tribbleStream, source);
        } else if (protocol.equals("ftp")) {
            org.broad.tribble.util.SeekableStream tribbleStream = new SeekableFTPStream(url);
            String source = url.toExternalForm();
            is = new SeekablePicardStream(tribbleStream, source);
        } else {
            throw new RuntimeException("Unknown protocol: " + protocol);
        }
        return is;
    }

    /**
     * Delete temporary files which are older than timeLimit.
     *
     * @param timeLimit Minimum age (in milliseconds) to delete. If null, default is 1 day
     * @throws IOException
     */
    public static void cleanTempDir(Long timeLimit) {
        if (timeLimit == null) {
            timeLimit = oneDay;
        }
        File dir = DirectoryManager.getCacheDirectory();
        File[] files = dir.listFiles();

        long time = System.currentTimeMillis();
        for (File f : files) {
            long age = time - f.lastModified();
            if (age > timeLimit) {
                f.delete();
            }
        }
    }


    // TODO -- revisit caching scheme,  do something for ftp loads
    File getIndexFile(URL url, String indexPath) throws IOException {

        String urlString = url.toString();
        indexFile = getTmpIndexFile(urlString);

        // Crude staleness check -- if more than a day old discard
        long age = System.currentTimeMillis() - indexFile.lastModified();
        if (age > oneDay) {
            indexFile.delete();
        }

        if (!indexFile.exists() || indexFile.length() < 1) {
            loadIndexFile(urlString, indexPath, indexFile);
            indexFile.deleteOnExit();
        }

        return indexFile;

    }

    private File getTmpIndexFile(String bamURL) throws IOException {
        File indexFile = indexFileCache.get(bamURL);
        if (indexFile == null) {
            indexFile = File.createTempFile("index_", ".bai", DirectoryManager.getCacheDirectory());
            indexFile.deleteOnExit();
            indexFileCache.put(bamURL, indexFile);
        }
        return indexFile;
    }

    private void loadIndexFile(String path, String indexPath, File indexFile) throws IOException {
        InputStream is = null;
        OutputStream os = null;

        try {
            String idx = (indexPath != null && indexPath.length() > 0) ? indexPath : path + ".bai";
            URL indexURL = new URL(idx);
            os = new FileOutputStream(indexFile);
            try {
                is = HttpUtils.getInstance().openConnectionStream(indexURL);
            } catch (FileNotFoundException e) {
                // Try other index convention
                String baseName = path.substring(0, path.length() - 4);
                indexURL = new URL(baseName + ".bai");

                try {
                    is = org.broad.igv.util.HttpUtils.getInstance().openConnectionStream(indexURL);
                } catch (FileNotFoundException e1) {
                    MessageUtils.showMessage("Index file not found for file: " + path);
                    throw new DataLoadException("Index file not found for file: " + path, path);
                }
            }
            byte[] buf = new byte[512000];
            int bytesRead;
            while ((bytesRead = is.read(buf)) != -1) {
                os.write(buf, 0, bytesRead);
            }

        } finally {
            if (is != null) {
                try {
                    is.close();
                } catch (IOException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }
            if (os != null) {
                try {
                    os.close();
                } catch (IOException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }

        }
    }


}
