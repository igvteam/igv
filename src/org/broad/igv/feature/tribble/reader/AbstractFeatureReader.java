/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.feature.tribble.reader;

import org.broad.igv.feature.tribble.TribbleIndexNotFoundException;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.*;
import org.broad.tribble.util.ParsingUtils;

import java.io.IOException;
import java.net.URL;
import java.util.Iterator;


/**
 * jrobinso
 * <p/>
 * the feature reader class, which uses indices and codecs to read in Tribble file formats.
 */
public abstract class AbstractFeatureReader<T extends Feature, SOURCE> implements org.broad.tribble.FeatureReader<T> {
    // the logging destination for this source
    //private final static Logger log = Logger.getLogger("BasicFeatureSource");

    // the path to underlying data source
    String path;

    // the query source, codec, and header
    // protected final QuerySource querySource;
    protected final FeatureCodec<T, SOURCE> codec;
    protected FeatureCodecHeader header;

    /**
     * Convenience overload which defaults to requiring an index.
     */
    public static <FEATURE extends Feature, SOURCE> AbstractFeatureReader<FEATURE, SOURCE>
    getFeatureReader(final ResourceLocator locator, final FeatureCodec<FEATURE, SOURCE> codec) throws TribbleException, TribbleIndexNotFoundException {

        try {
            // Test for tabix index
            if (isTabix(locator)) {
                if (!(codec instanceof AsciiFeatureCodec)) {
                    throw new TribbleException("Tabix indexed files only work with ASCII codecs, but received non-Ascii codec " +
                            codec.getClass().getSimpleName());
                }
                return new TabixFeatureReader<FEATURE, SOURCE>(locator, (AsciiFeatureCodec) codec);
            }
            // Not tabix => tribble index file (might be gzipped, but not block gzipped)
            else {
                return new TribbleFeatureReader<FEATURE, SOURCE>(locator, codec);
            }
        } catch (IOException e) {
            throw new TribbleException.MalformedFeatureFile("Unable to create BasicFeatureReader using feature file ", locator.getPath(), e);
        } catch (TribbleException e) {
            e.setSource(locator.getPath());
            throw e;
        }
    }

    private static boolean isTabix(ResourceLocator locator) throws IOException {
        String tabxIndex = locator.getIndexPath();
        if (tabxIndex == null) {
            if (HttpUtils.isRemoteURL(locator.getPath())) {
                final URL url = new URL(locator.getPath());
                String path = url.getPath();
                String indexPath = path + ".tbi";   // Strip off parameters
                tabxIndex = locator.getPath().replace(path, indexPath);
            } else {
                tabxIndex = locator.getPath() + ".tbi";
            }
        }
        boolean isTabix =  locator.getPath().endsWith(".gz") && FileUtils.resourceExists(tabxIndex);
        if(isTabix) {
            locator.setIndexPath(tabxIndex);
        }
        return isTabix;
    }

    protected AbstractFeatureReader(final String path, final FeatureCodec<T, SOURCE> codec) {
        this.path = path;
        this.codec = codec;
    }


    /**
     * get the header
     *
     * @return the header object we've read-in
     */
    public Object getHeader() {
        return header.getHeaderValue();
    }

    abstract public boolean hasIndex();

    static class EmptyIterator<T extends Feature> implements CloseableTribbleIterator<T> {
        public Iterator iterator() {
            return this;
        }

        public boolean hasNext() {
            return false;
        }

        public T next() {
            return null;
        }

        public void remove() {
        }

        @Override
        public void close() {
        }
    }
}
