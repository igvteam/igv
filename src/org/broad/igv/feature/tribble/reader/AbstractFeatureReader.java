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

import org.broad.tribble.*;
import org.broad.tribble.index.Index;
import org.broad.tribble.util.ParsingUtils;

import java.io.IOException;
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

    /** Convenience overload which defaults to requiring an index. */
    public static <FEATURE extends Feature, SOURCE> AbstractFeatureReader<FEATURE, SOURCE> getFeatureReader(final String featureFile, final FeatureCodec<FEATURE, SOURCE> codec) throws TribbleException {
        return getFeatureReader(featureFile, codec, true);
    }

    /**
     * factory for unknown file type,  could be ascii, binary, or could be tabix, or something else.
     *
     * @param featureResource the feature file to create from
     * @param codec           the codec to read features with
     */
    public static <FEATURE extends Feature, SOURCE> AbstractFeatureReader<FEATURE, SOURCE> getFeatureReader(final String featureResource, final FeatureCodec<FEATURE, SOURCE> codec, final boolean requireIndex) throws TribbleException {

        try {
            // Test for tabix index
            if (featureResource.endsWith(".gz") && ParsingUtils.resourceExists(featureResource + ".tbi")) {
                if ( ! (codec instanceof AsciiFeatureCodec) )
                    throw new TribbleException("Tabix indexed files only work with ASCII codecs, but received non-Ascii codec " + codec.getClass().getSimpleName());
                return new TabixFeatureReader<FEATURE, SOURCE>(featureResource, (AsciiFeatureCodec) codec);
            }
            // Not tabix => tribble index file (might be gzipped, but not block gzipped)
            else {
                return new TribbleIndexedFeatureReader<FEATURE, SOURCE>(featureResource, codec, requireIndex);
            }
        } catch (IOException e) {
            throw new TribbleException.MalformedFeatureFile("Unable to create BasicFeatureReader using feature file ", featureResource, e);
        } catch (TribbleException e) {
            e.setSource(featureResource);
            throw e;
        }
    }

    /**
     * Return a reader with a supplied index.
     *
     * @param featureResource the path to the source file containing the features
     * @param codec used to decode the features
     * @param index index of featureResource
     * @return a reader for this data
     * @throws TribbleException
     */
    public static <FEATURE extends Feature, SOURCE> AbstractFeatureReader<FEATURE, SOURCE> getFeatureReader(final String featureResource, final FeatureCodec<FEATURE, SOURCE>  codec, final Index index) throws TribbleException {
        try {
            return new TribbleIndexedFeatureReader<FEATURE, SOURCE>(featureResource, codec, index);
        } catch (IOException e) {
            throw new TribbleException.MalformedFeatureFile("Unable to create AbstractFeatureReader using feature file ", featureResource, e);
        }

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

    static class EmptyIterator<T extends Feature> implements CloseableTribbleIterator<T> {
        public Iterator iterator() { return this; }
        public boolean hasNext() { return false; }
        public T next() { return null; }
        public void remove() { }
        @Override public void close() { }
    }
}
