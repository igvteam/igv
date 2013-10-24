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

package org.broad.igv.track;

import org.broad.igv.data.DataSource;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.*;
import org.broad.igv.feature.tribble.reader.AbstractFeatureReader;
import org.broad.igv.tdf.TDFDataSource;
import org.broad.igv.tdf.TDFReader;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.RuntimeUtils;
import org.broad.tribble.CloseableTribbleIterator;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.FeatureReader;

import java.io.IOException;
import java.util.*;

/**
 * @author jrobinso
 * @date Jun 27, 2010
 */
public class TribbleFeatureSource implements org.broad.igv.track.FeatureSource {

    IGVFeatureReader reader;
    DataSource coverageSource;
    boolean isVCF;
    Genome genome;
    boolean useCache;

    /**
     * Map of IGV chromosome name -> source name
     */
    Map<String, String> chrNameMap = new HashMap<String, String>();
    private int featureWindowSize;
    Object header;
    Class featureClass;

    public TribbleFeatureSource(String path, Genome genome) throws IOException {
        this(path, genome, true);
    }


    public TribbleFeatureSource(String path, Genome genome, boolean useCache) throws IOException {
        this.genome = genome;
        this.useCache = useCache;
        init(path);
    }

    private void initCoverageSource(String covPath) {
        if (ParsingUtils.pathExists(covPath)) {
            TDFReader reader = TDFReader.getReader(covPath);
            coverageSource = new TDFDataSource(reader, 0, "", genome);
        }
    }

    protected void init(String path) {

        FeatureCodec codec = CodecFactory.getCodec(path, genome);
        isVCF = codec.getClass() == VCFWrapperCodec.class;
        featureClass = codec.getFeatureType();
        AbstractFeatureReader basicReader = AbstractFeatureReader.getFeatureReader(path, codec, true);
        header = basicReader.getHeader();
        initFeatureWindowSize(basicReader);
        reader = useCache ?
                new CachingFeatureReader(basicReader, 5, getFeatureWindowSize()) :
                new TribbleReaderWrapper(basicReader);

        if (genome != null) {
            Collection<String> seqNames = reader.getSequenceNames();
            if (seqNames != null)
                for (String seqName : seqNames) {
                    String igvChr = genome.getChromosomeAlias(seqName);
                    if (igvChr != null && !igvChr.equals(seqName)) {
                        chrNameMap.put(igvChr, seqName);
                    }
                }
        }

        initCoverageSource(path + ".tdf");
    }

    public Class getFeatureClass() {
        return featureClass;
    }

    public Iterator<Feature> getFeatures(String chr, int start, int end) throws IOException {

        String seqName = chrNameMap.get(chr);
        if (seqName == null) seqName = chr;

        return reader.query(seqName, start, end);
    }

    /**
     * Return coverage values overlapping the query interval.   At this time Tribble sources do not provide
     * coverage values
     *
     * @param chr
     * @param start
     * @param end
     * @param zoom
     * @return
     */
    public List<LocusScore> getCoverageScores(String chr, int start, int end, int zoom) {
        return coverageSource == null ? null :
                coverageSource.getSummaryScoresForRange(chr, start, end, zoom);
    }

    public int getFeatureWindowSize() {
        return featureWindowSize;
    }

    public void setFeatureWindowSize(int size) {
        this.featureWindowSize = size;
        if (reader instanceof CachingFeatureReader) {
            ((CachingFeatureReader) reader).setBinSize(size);
        }
    }

    public Object getHeader() {
        return header;
    }


    /**
     * Estimate an appropriate feature window size.
     *
     * @param reader
     */
    private void initFeatureWindowSize(FeatureReader reader) {

        CloseableTribbleIterator<org.broad.tribble.Feature> iter = null;

        try {
            double mem = RuntimeUtils.getAvailableMemory();
            iter = reader.iterator();
            if (iter.hasNext()) {

                int nSamples = isVCF ? 100 : 1000;
                org.broad.tribble.Feature firstFeature = iter.next();
                org.broad.tribble.Feature lastFeature = firstFeature;
                String chr = firstFeature.getChr();
                int n = 1;
                long len = 0;
                while (iter.hasNext() && n < nSamples) {
                    org.broad.tribble.Feature f = iter.next();
                    if (f != null) {
                        n++;
                        if (f.getChr().equals(chr)) {
                            lastFeature = f;
                        } else {
                            len += lastFeature.getEnd() - firstFeature.getStart() + 1;
                            firstFeature = f;
                            lastFeature = f;
                            chr = f.getChr();
                        }

                    }
                }
                double dMem = mem - RuntimeUtils.getAvailableMemory();
                double bytesPerFeature = Math.max(100, dMem / n);

                len += lastFeature.getEnd() - firstFeature.getStart() + 1;
                double featuresPerBase = ((double) n) / len;

                double targetBinMemory = 20000000;  // 20  mega bytes
                int maxBinSize = isVCF ? 1000000 : Integer.MAX_VALUE;
                int bs = Math.min(maxBinSize, (int) (targetBinMemory / (bytesPerFeature * featuresPerBase)));
                featureWindowSize = Math.max(1000000, bs);
            } else {
                featureWindowSize = Integer.MAX_VALUE;
            }
        } catch (IOException e) {
            featureWindowSize = 1000000;
        } finally {
            if (iter != null) iter.close();
        }
    }

    public Collection<String> getChrNames() {
        return chrNameMap.keySet();
    }
}
