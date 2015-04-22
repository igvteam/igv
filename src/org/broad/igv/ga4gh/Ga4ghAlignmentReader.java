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

    public Ga4ghAlignmentReader(Ga4ghProvider provider, String readsetId) {
        this.provider = provider;
        this.readsetId = readsetId;
    }

    @Override
    public void close() throws IOException {

    }

    @Override
    public List<String> getSequenceNames() {

        if (sequenceNames == null) {
            try {
                loadMetadata();
            } catch (Exception e) {
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

    @Override
    public CloseableIterator<Alignment> query(String sequence, int start, int end, boolean contained) throws IOException {

        List<Alignment> alignmentList = Ga4ghAPIHelper.searchReads(provider, readsetId, sequence, start, end);

        return alignmentList == null ? null : new MIterator(alignmentList);
    }

    @Override
    public boolean hasIndex() {
        return true;
    }

    private void loadMetadata() throws IOException {

        String authKey = provider.authKey;
        String baseURL = provider.baseURL;
        URL url = new URL(baseURL + "/readgroupsets/" + readsetId + (authKey == null ? "" : "?key=" + authKey));   // TODO -- field selection?

        String result = HttpUtils.getInstance().getContentsAsString(url);
        JsonParser parser = new JsonParser();
        JsonObject root = parser.parse(result).getAsJsonObject();

        if (root.has("referenceSetId")) {
            String referenceSetId = root.getAsJsonPrimitive("referenceSetId").getAsString();

            List<JsonObject> refererences = Ga4ghAPIHelper.searchReferences(provider, referenceSetId, 1000);

            sequenceNames = new ArrayList();

            for (JsonObject refObject : refererences) {
                sequenceNames.add(refObject.getAsJsonPrimitive("name").getAsString());
            }
        }
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

}
