package org.broad.igv.bigwig;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import org.apache.log4j.Logger;
import org.broad.igv.data.BasicScore;
import org.broad.igv.data.DataSource;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackType;
import org.broad.igv.track.WindowFunction;
import org.broad.igv.util.HttpUtils;

import java.io.IOException;
import java.net.URL;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.*;

import static org.broad.igv.data.AbstractDataSource.ORDERED_WINDOW_FUNCTIONS;

public class D4ServerSource implements DataSource {

    private static Logger log = Logger.getLogger(DataSource.class);

    String url;
    double min = 0;
    double max = 0.001;
    Map<String, String> chrAliasTable;
    Genome genome;

    public D4ServerSource(String url, Genome genome) {
        this.url = url;
        this.genome = genome;
    }

    private Map<String, String> getChrAliasTable() {
        if (chrAliasTable == null && genome != null) {
            chrAliasTable = new HashMap<>();
            try {
                String queryURL = this.url.replace("d4get://", "http://") + "?class=header";
                String json = HttpUtils.getInstance().getContentsAsJSON(new URL(queryURL));
                JsonParser parser = new JsonParser();
                JsonArray obj = parser.parse(json).getAsJsonArray();
                Iterator<JsonElement> iter = obj.iterator();
                while (iter.hasNext()) {
                    String chrName = iter.next().getAsString();
                    String canonicalName = genome.getCanonicalChrName(chrName);
                    chrAliasTable.put(canonicalName, chrName);
                }
            } catch (IOException e) {
                log.error("Error initializing chromosomes", e);
            }
        }
        return chrAliasTable;
    }

    @Override
    public double getDataMax() {
        return min;
    }

    @Override
    public double getDataMin() {
        return max;
    }

    @Override
    public List<LocusScore> getSummaryScoresForRange(String chr, int start, int end, int zoom) {
        try {

            Map<String, String> chrAliasTable = getChrAliasTable();
            String queryChr = chrAliasTable != null && chrAliasTable.containsKey(chr) ?
                    chrAliasTable.get(chr) : chr;


            String queryURL = this.url.replace("d4get://", "http://") + "?chr=" + queryChr + "&start=" + start + "&end=" + end;
            byte[] bytes = HttpUtils.getInstance().getContentsAsBytes(new URL(queryURL), null);
            ByteBuffer byteBuffer = ByteBuffer.wrap(bytes);
            byteBuffer.order(ByteOrder.LITTLE_ENDIAN);
            int dataStart = byteBuffer.getInt();
            int stepSize = byteBuffer.getInt();
            int nPoints = byteBuffer.getInt();
            List<LocusScore> scores = new ArrayList<>(nPoints);
            for (int i = 0; i < nPoints; i++) {
                float value = byteBuffer.getFloat();
                scores.add(new BasicScore(dataStart, dataStart + stepSize, value));
                dataStart += stepSize;
                min = Math.min(min, value);
                max = Math.max(max, value);

            }
            return scores;
        } catch (IOException e) {
            log.error("Error reading from D4 server", e);
            return null;
        }
    }

    @Override
    public TrackType getTrackType() {
        return null;
    }

    @Override
    public void setWindowFunction(WindowFunction statType) {

    }

    @Override
    public boolean isLogNormalized() {
        return false;
    }

    @Override
    public WindowFunction getWindowFunction() {
        return null;
    }

//    @Override
//    public Collection<WindowFunction> getAvailableWindowFunctions() {
//        return ORDERED_WINDOW_FUNCTIONS;
//    }

    @Override
    public void dispose() {

    }
}
