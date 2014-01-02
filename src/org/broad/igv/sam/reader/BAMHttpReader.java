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
import net.sf.samtools.seekablestream.SeekableStream;
import net.sf.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.sam.SamAlignment;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.stream.IGVSeekableBufferedStream;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

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
public class BAMHttpReader implements AlignmentReader<SamAlignment> {

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
            indexFile = getIndexFile(locator);
            if (indexFile == null) {
                throw new RuntimeException("Could not load index file for file: " + url.getPath());
            }

            SeekableStream ss = new IGVSeekableBufferedStream(IGVSeekableStreamFactory.getInstance().getStreamFor(url), 128000);
            //SeekableStream ss = getSeekableStream(url);
            log.debug("Initializing SAMFileReader");

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

    public SAMFileHeader getFileHeader() {
        if (header == null) {
            header = reader.getFileHeader();
        }
        return header;
    }

    public Set<String> getPlatforms() {
        return AlignmentReaderFactory.getPlatforms(getFileHeader());
    }

    public boolean hasIndex() {
        return indexFile != null && indexFile.exists();
    }

    public List<String> getSequenceNames() {
        if (sequenceNames == null) {
            SAMFileHeader header = getFileHeader();
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


    public CloseableIterator<SamAlignment> iterator() {
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

    public CloseableIterator<SamAlignment> query(String sequence, int start, int end, boolean contained) {
        try {
            if (reader == null) {
                SeekableStream ss = new IGVSeekableBufferedStream(IGVSeekableStreamFactory.getInstance().getStreamFor(url));
                reader = new SAMFileReader(ss, indexFile, false);
            }
            CloseableIterator<SAMRecord> iter = reader.query(sequence, start + 1, end, contained);
            return new WrappedIterator(iter);
        } catch (IOException e) {
            log.error("Error opening SAM reader", e);
            throw new RuntimeException("Error opening SAM reader", e);
        }
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

    File getIndexFile(ResourceLocator locator) throws IOException {

        log.debug("Getting index for " + url + ". Index path " + locator.getBamIndexPath());
        String urlString = url.toString();
        indexFile = getTmpIndexFile(urlString);

        // Crude staleness check -- if more than a day old discard
        long age = System.currentTimeMillis() - indexFile.lastModified();
        if (age > oneDay) {
            indexFile.delete();
        }

        if (!indexFile.exists() || indexFile.length() < 1) {
            loadIndexFile(locator.getBamIndexPath(), indexFile);
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

    private void loadIndexFile(String indexPath, File indexFile) throws IOException {
        InputStream is = null;
        OutputStream os = null;

        try {
            URL indexURL = new URL(indexPath);
            os = new FileOutputStream(indexFile);
            boolean foundIndex = true;
            try {
                is = HttpUtils.getInstance().openConnectionStream(indexURL);
            } catch (FileNotFoundException e) {
                // Try other index convention
                indexPath = indexPath.replace(".bam.bai", ".bai");
                indexURL = new URL(indexPath);
                try {
                    is = org.broad.igv.util.HttpUtils.getInstance().openConnectionStream(indexURL);
                } catch (FileNotFoundException e1) {

                    if (!Globals.isHeadless() && IGV.hasInstance()) {
                        String tmp = MessageUtils.showInputDialog("Index file not found. Enter path to index file", indexPath);
                        if (tmp != null) {
                            try {
                                indexURL = new URL(tmp);
                                is = org.broad.igv.util.HttpUtils.getInstance().openConnectionStream(indexURL);
                            } catch (FileNotFoundException e2) {
                                foundIndex = false;
                            }
                        } else {
                            foundIndex = false;
                        }
                    }
                }

            }
            if (!foundIndex) {
                String msg = "Index file not found: " + indexPath;
                throw new DataLoadException(msg, indexPath);
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
                    log.error(e.getMessage(), e);
                }
            }
            if (os != null) {
                try {
                    os.close();
                } catch (IOException e) {
                    log.error(e.getMessage(), e);
                }
            }

        }
    }


}
