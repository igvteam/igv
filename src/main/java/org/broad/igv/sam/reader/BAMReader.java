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

package org.broad.igv.sam.reader;

import htsjdk.samtools.*;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.sam.EmptyAlignmentIterator;
import org.broad.igv.sam.PicardAlignment;
import org.broad.igv.sam.cram.IGVReferenceSource;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.URLUtils;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Sep 22, 2009
 * Time: 2:21:04 PM
 */
public class BAMReader implements AlignmentReader<PicardAlignment> {

    static Logger log = Logger.getLogger(BAMReader.class);

    private final ResourceLocator locator;

    SAMFileHeader header;
    htsjdk.samtools.SamReader reader;
    List<String> sequenceNames;
    private boolean indexed = false; // False until proven otherwise
    private Map<String, Long> sequenceDictionary;

    public BAMReader(ResourceLocator locator, boolean requireIndex) throws IOException {
        this.locator = locator;
        reader = getSamReader(locator, requireIndex);
        header = reader.getFileHeader();
    }

    private SamReader getSamReader(ResourceLocator locator, boolean requireIndex) throws IOException {

        boolean isLocal = locator.isLocal();
        final SamReaderFactory factory = SamReaderFactory.makeDefault().
                referenceSource(new IGVReferenceSource()).
                validationStringency(ValidationStringency.SILENT);
        SamInputResource resource;

        if (isLocal) {
            resource = SamInputResource.of(new File(locator.getPath()));
        } else {
            URL url = HttpUtils.createURL(locator.getPath());
            if (requireIndex) {
                resource = SamInputResource.of(IGVSeekableStreamFactory.getInstance().getStreamFor(url));
            } else {
                resource = SamInputResource.of(new BufferedInputStream(HttpUtils.getInstance().openConnectionStream(url)));
            }
        }

        if (requireIndex) {

            String indexPath = getExplicitIndexPath(locator);
            if (indexPath == null || indexPath.length() == 0) {
                indexPath = getIndexPath(locator.getPath());
            }

            indexed = true;
            if (isLocal) {
                File indexFile = new File(indexPath);
                resource = resource.index(indexFile);
            } else {
                SeekableStream indexStream = IGVSeekableStreamFactory.getInstance().getStreamFor(HttpUtils.createURL(indexPath));
                resource = resource.index(indexStream);
            }
        }

        return factory.open(resource);

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

    public boolean hasIndex() {
        return indexed;
    }

    public Set<String> getPlatforms() {
        return AlignmentReaderFactory.getPlatforms(getFileHeader());
    }


    public List<String> getSequenceNames() {
        if (sequenceNames == null) {
            loadSequenceDictionary();
        }
        return sequenceNames;
    }

    @Override
    public Map<String, Long> getSequenceDictionary() {
        if (sequenceDictionary == null) {
            loadSequenceDictionary();
        }
        return sequenceDictionary;
    }

    private void loadSequenceDictionary() {

        SAMFileHeader header = getFileHeader();
        sequenceNames = new ArrayList();
        sequenceDictionary = new HashMap<>();

        List<SAMSequenceRecord> records = header.getSequenceDictionary().getSequences();
        if (records.size() > 0) {
            for (SAMSequenceRecord rec : header.getSequenceDictionary().getSequences()) {
                String chr = rec.getSequenceName();
                Long size = new Long(rec.getSequenceLength());

                sequenceNames.add(chr);
                sequenceDictionary.put(chr, size);
            }
        }
    }


    public CloseableIterator<PicardAlignment> iterator() {
        return new WrappedIterator(reader.iterator());
    }

    public CloseableIterator<PicardAlignment> query(String sequence, int start, int end, boolean contained) {
        if (sequenceDictionary != null && !sequenceDictionary.containsKey(sequence)) {
            return EMPTY_ITERATOR;
        } else {
            CloseableIterator<SAMRecord> iter = null;
            try {
                iter = reader.query(sequence, start + 1, end, contained);
            } catch (IllegalArgumentException e) {
                log.error("Error querying for sequence: " + sequence, e);
                return new EmptyAlignmentIterator();
            }
            return new ListIterator(iter);
        }
    }

    /**
     * Fetch an explicitly set index path, either via the ResourceLocator or as a parameter in a URL
     *
     * @param locator
     * @return the index path, or null if no index path is set
     */
    private String getExplicitIndexPath(ResourceLocator locator) {

        String p = locator.getPath().toLowerCase();
        String idx = locator.getIndexPath();

        if (idx == null && (p.startsWith("http://") || p.startsWith("https://"))) {
            try {
                URL url = HttpUtils.createURL(locator.getPath());
                String queryString = url.getQuery();
                if (queryString != null) {
                    Map<String, String> parameters = URLUtils.parseQueryString(queryString);
                    if (parameters.containsKey("index")) {
                        idx = parameters.get("index");

                    }
                }
            } catch (MalformedURLException e) {
                log.error("Error parsing url: " + locator.getPath());
            }
        }
        return idx;
    }

    /**
     * Try to guess the index path.
     *
     * @param pathOrURL
     * @return
     * @throws IOException
     */
    private String getIndexPath(String pathOrURL) throws IOException {

        List<String> pathsTried = new ArrayList<String>();

        String indexPath;

        if (URLUtils.isURL(pathOrURL)) {

            String path = URLUtils.getPath(pathOrURL);

            if (path.endsWith(".bam")) {
                // Try .bam.bai
                indexPath = URLUtils.addExtension(pathOrURL, ".bai");
                pathsTried.add(indexPath);
                if (HttpUtils.getInstance().resourceAvailable(indexPath)) {
                    return indexPath;
                }

                // Try .bai
                indexPath = URLUtils.replaceExtension(pathOrURL, ".bam", ".bai");
                pathsTried.add(indexPath);
                if (HttpUtils.getInstance().resourceAvailable(indexPath)) {
                    return indexPath;
                }

                // Try .bam.csi
                indexPath = URLUtils.addExtension(pathOrURL, ".csi");
                pathsTried.add(indexPath);
                if (HttpUtils.getInstance().resourceAvailable(indexPath)) {
                    return indexPath;
                }

                // Try .csi
                indexPath = URLUtils.replaceExtension(pathOrURL, ".bam", ".csi");
                pathsTried.add(indexPath);
                if (HttpUtils.getInstance().resourceAvailable(indexPath)) {
                    return indexPath;
                }
            }

            //  cram
            if (path.endsWith(".cram")) {
                indexPath = URLUtils.addExtension(pathOrURL, ".crai");
                if (FileUtils.resourceExists(indexPath)) {
                    pathsTried.add(indexPath);
                    return indexPath;
                } else {
                    indexPath = pathOrURL.substring(0, pathOrURL.length() - 5) + ".crai";
                    if (FileUtils.resourceExists(indexPath)) {
                        return indexPath;
                    }
                }
            }


        } else {
            // Local file

            indexPath = pathOrURL + ".bai";

            if (FileUtils.resourceExists(indexPath)) {
                return indexPath;
            }

            if (indexPath.contains(".bam.bai")) {
                indexPath = indexPath.replaceFirst(".bam.bai", ".bai");
                pathsTried.add(indexPath);
                if (FileUtils.resourceExists(indexPath)) {
                    return indexPath;
                }
            } else {
                indexPath = indexPath.replaceFirst(".bai", ".bam.bai");
                pathsTried.add(indexPath);
                if (FileUtils.resourceExists(indexPath)) {
                    return indexPath;
                }
            }


            // Try .bam.csi
            indexPath = pathOrURL + ".csi";
            pathsTried.add(indexPath);
            if (FileUtils.resourceExists(indexPath)) {
                return indexPath;
            }

            // Try .csi
            if (pathOrURL.endsWith(".bam")) {
                indexPath = pathOrURL.substring(0, pathOrURL.length() - 4) + ".csi";
                pathsTried.add(indexPath);
                if (FileUtils.resourceExists(indexPath)) {
                    return indexPath;
                }
            }

            if (pathOrURL.endsWith(".cram")) {
                indexPath = pathOrURL + ".crai";
                if (FileUtils.resourceExists(indexPath)) {
                    return indexPath;
                } else {
                    indexPath = pathOrURL.substring(0, pathOrURL.length() - 5) + ".crai";
                    if (FileUtils.resourceExists(indexPath)) {
                        return indexPath;
                    }
                }
            }


        }

        String defaultValue = pathOrURL + (pathOrURL.endsWith(".cram") ? ".crai" : ".bai");
        indexPath = MessageUtils.showInputDialog(
                "Index is required, but no index found.  Please enter path to index file:",
                defaultValue);
        if (indexPath != null && FileUtils.resourceExists(indexPath)) {
            return indexPath;
        }


        String msg = "Index file not found.  Tried ";
        for (String p : pathsTried) {
            msg += "<br>" + p;
        }
        throw new DataLoadException(msg, indexPath);

    }


    static CloseableIterator<PicardAlignment> EMPTY_ITERATOR = new CloseableIterator<PicardAlignment>() {
        @Override
        public void close() {

        }

        @Override
        public boolean hasNext() {
            return false;
        }

        @Override
        public PicardAlignment next() {
            return null;
        }
    };

    static class ListIterator implements CloseableIterator<PicardAlignment> {

        Iterator<PicardAlignment> iterator;

        public ListIterator(CloseableIterator<SAMRecord> iter) {
            List<PicardAlignment> alignments = new ArrayList<>(1000);
            while (iter.hasNext()) {
                alignments.add(new PicardAlignment(iter.next()));
            }
            iter.close();
            iterator = alignments.iterator();

        }

        public void close() {
            iterator = null;
        }

        public boolean hasNext() {
            return iterator.hasNext();
        }

        public PicardAlignment next() {
            return iterator.next();
        }

        public void remove() {
            iterator = null;
        }


    }


}
