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

package org.broad.igv.ga4gh;

import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.util.HttpUtils;

import java.io.IOException;
import java.net.URL;
import java.util.*;

/**
 * Class for testing GlobalAlliance API.  Loads json from a text file, for development only
 * <p/>
 * Created by jrobinso on 7/18/14.
 */

public class Ga4ghAlignmentReader implements AlignmentReader<Alignment> {

    private static Logger log = Logger.getLogger(Ga4ghAlignmentReader.class);

    String readsetId;
    List<String> sequenceNames;
    Ga4ghProvider provider;
    Map<String, String> chromosomeMappings;

    public Ga4ghAlignmentReader(Ga4ghProvider provider, String readsetId) {
        this.provider = provider;
        this.readsetId = readsetId;
        this.chromosomeMappings = new HashMap<String, String>();
    }

    @Override
    public void close() throws IOException {

    }

    @Override
    public List<String> getSequenceNames() throws IOException {

        if (sequenceNames == null) {
            try {
                loadMetadata();
            }
            catch (IOException e) {
                throw e;
            }
            catch (Exception e) {
                log.error("Error fetching metadata", e);
            }
        }
        return sequenceNames;
    }

    @Override
    public SAMFileHeader getFileHeader() {
        return null;
    }

    @Override
    public Set<String> getPlatforms() {
        return null;
    }

    @Override
    public CloseableIterator<Alignment> iterator() {
        throw new RuntimeException("Iterating over ga4gh datasets is not supported");
    }

    /**
     * Return an iterator over alignments overlapping the given region.  We first try the sequence name as given,
     * then try an alias if we have one.  This is to handle the common 1 vs chr1 naming problem stemming from
     * competing conventions.
     *
     * @param sequence
     * @param start     0-based start location
     * @param end       0-based, exclusive-end coordinate
     * @param contained
     * @return
     * @throws IOException
     */
    @Override
    public CloseableIterator<Alignment> query(String sequence, int start, int end, boolean contained) throws IOException {

        boolean sequenceNameConfirmed = false;
        if (chromosomeMappings.containsKey(sequence)) {
            sequence = chromosomeMappings.get(sequence);
            sequenceNameConfirmed = true;
        }

        boolean handleError = sequenceNameConfirmed;

        List<Alignment> alignmentList = Ga4ghAPIHelper.searchReads(provider, readsetId, sequence, start, end, handleError);

        if (alignmentList != null) {
            chromosomeMappings.put(sequence, sequence);
            return new MIterator(alignmentList);
        } else if (sequenceNameConfirmed == false && COMMON_ALIASES.containsKey(sequence)) {
            String altSequence = COMMON_ALIASES.get(sequence);
            log.info("Sequence \"" + sequence + "\" failed.  Trying \"" + altSequence + "\"");
            alignmentList = Ga4ghAPIHelper.searchReads(provider, readsetId, altSequence, start, end, true);
            if (alignmentList != null) {
                chromosomeMappings.put(sequence, altSequence);
            }
            return alignmentList == null ? null : new MIterator(alignmentList);
        } else {
            return null;
        }

    }

    @Override
    public boolean hasIndex() {
        return true;
    }

    private void loadMetadata() throws IOException {

        String authKey = provider.apiKey;
        String baseURL = provider.baseURL;
        URL url = new URL(baseURL + "/readgroupsets/" + readsetId + (authKey == null ? "" : "?key=" + authKey));   // TODO -- field selection?

        String result = HttpUtils.getInstance().getContentsAsString(url);
        JsonParser parser = new JsonParser();
        JsonObject root = parser.parse(result).getAsJsonObject();

//        if (root.has("referenceSetId")) {
//            String referenceSetId = root.getAsJsonPrimitive("referenceSetId").getAsString();
//
//            List<JsonObject> refererences = Ga4ghAPIHelper.searchReferences(provider, referenceSetId, 1000);
//
//            sequenceNames = new ArrayList();
//
//            for (JsonObject refObject : refererences) {
//                sequenceNames.add(refObject.getAsJsonPrimitive("name").getAsString());
//            }
//        }
    }

    public static boolean supportsFileType(String type) {
        return type.equals(Ga4ghAPIHelper.RESOURCE_TYPE);
    }

    class MIterator implements CloseableIterator<Alignment> {

        Iterator<Alignment> iter;

        MIterator(List<Alignment> alignmentList) {
            iter = alignmentList.iterator();
        }

        @Override
        public void close() {

        }

        @Override
        public boolean hasNext() {
            return iter.hasNext();
        }

        @Override
        public Alignment next() {
            return iter.next();
        }

        @Override
        public void remove() {
            iter.remove();
        }
    }


    static Map<String, String> COMMON_ALIASES = new HashMap<String, String>();

    static {
        for (int i = 1; i <= 23; i++) {
            COMMON_ALIASES.put(("" + i), ("chr" + i));
            COMMON_ALIASES.put(("chr" + i), ("" + i));
        }
        COMMON_ALIASES.put("chrX", "X");
        COMMON_ALIASES.put("X", "chrX");
        COMMON_ALIASES.put("chrY", "Y");
        COMMON_ALIASES.put("Y", "chrY");
        COMMON_ALIASES.put("chrM", "MT");
        COMMON_ALIASES.put("MT", "chrM");
    }

}
