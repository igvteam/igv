package org.broad.igv.sam.uncalled4;

import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import org.broad.igv.sam.reader.AlignmentReader;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Uncalled4Utils {

    public static Map<String, Map<String, Float>> getScaleFactors(AlignmentReader reader) {

        List<String> comments = reader.getFileHeader().getComments();
        for (String c : comments) {
            if (c.startsWith("@CO\tUNC:")) {
                Map<String, Map<String, Float>> trackMap = new HashMap<>();
                String jsonString = c.substring(8);
                JsonElement element = JsonParser.parseString(jsonString);
                JsonObject tracks = element.getAsJsonObject().get("tracks").getAsJsonObject();
                for(String t : tracks.keySet()) {
                    JsonObject layers = tracks.get(t).getAsJsonObject().get("layers").getAsJsonObject();
                    float ucScale = layers.get("uc").getAsJsonObject().get("scale").getAsFloat();
                    float udScale = layers.get("ud").getAsJsonObject().get("scale").getAsFloat();
                    Map<String, Float> scales = new HashMap<>();
                    scales.put("uc", ucScale);
                    scales.put("ud", udScale);
                    trackMap.put(t, scales);
                }
                return trackMap;
            }
        }
        return null;
    }



}
