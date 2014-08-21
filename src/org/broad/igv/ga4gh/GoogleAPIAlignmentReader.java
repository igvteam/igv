package org.broad.igv.ga4gh;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.CloseableIterator;
import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
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

public class GoogleAPIAlignmentReader implements AlignmentReader<Alignment> {

    private static Logger log = Logger.getLogger(GoogleAPIAlignmentReader.class);

    String readsetId;
    List<String> sequenceNames;

    public GoogleAPIAlignmentReader(String readsetId) {
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
            } catch (IOException e) {
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

        List<Alignment> alignmentList =  GoogleAPIHelper.reads(readsetId, sequence, start, end);

        return new MIterator(alignmentList);
    }

    @Override
    public boolean hasIndex() {
        return true;
    }

    private void loadMetadata() throws IOException {

        String authKey = PreferenceManager.getInstance().get(PreferenceManager.GOOGLE_API_KEY);
        String baseURL = PreferenceManager.getInstance().get(PreferenceManager.GOOGLE_BASE_URL);
        URL url = new URL(baseURL + "/readsets/" + readsetId + "?key=" + authKey);   // TODO -- field selection?

        String result = HttpUtils.getInstance().getContentsAsString(url);
        JsonParser parser = new JsonParser();
        JsonObject root = parser.parse(result).getAsJsonObject();

        sequenceNames = new ArrayList();
        JsonArray fileData = root.getAsJsonArray("fileData");
        Iterator<JsonElement> fileIter = fileData.iterator();
        while (fileIter.hasNext()) {

            JsonObject fileObject = fileIter.next().getAsJsonObject();
            JsonArray refSequences = fileObject.getAsJsonArray("refSequences");

            Iterator<JsonElement> iter = refSequences.iterator();
            while (iter.hasNext()) {
                sequenceNames.add(iter.next().getAsJsonObject().get("name").getAsString());
            }
        }
    }

    public static boolean supportsFileType(String type) {
        return type.equals(GoogleAPIHelper.RESOURCE_TYPE);
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
