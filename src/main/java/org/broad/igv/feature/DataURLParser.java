package org.broad.igv.feature;

import htsjdk.tribble.Feature;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.util.ResourceLocator;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.net.URLDecoder;
import java.net.URLEncoder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

public class DataURLParser  {

    List<Feature> features;
    TrackProperties trackProperties;

    public DataURLParser() {
    }

    public  List<Feature> parseFeatures(String dataURL, String format, Genome genome) {

        if (!dataURL.startsWith("data:")) {
            throw new Error("Not a dataURL: " + dataURL);
        }

        int commaIndex = dataURL.indexOf(',');
        if (commaIndex < 0) {
            throw new Error("dataURL missing commas");
        }

        // TODO -- check optional media type - only text/plain supported

        String contents = URLDecoder.decode(dataURL.substring(commaIndex + 1), StandardCharsets.UTF_8);
        ResourceLocator locator = new ResourceLocator("");
        locator.setFormat(format);
        FeatureParser fp = AbstractFeatureParser.getInstanceFor(locator, genome);
        if (fp == null) {
            throw new Error("Unsupported format: " + format);
        }

        BufferedReader reader = new BufferedReader(new StringReader(contents));

        features = fp.loadFeatures(reader, genome);
        trackProperties = fp.getTrackProperties();

        return features;
    }

    public List<Feature> getFeatures() {
        return features;
    }

    public TrackProperties getTrackProperties() {
        return trackProperties;
    }

    /**
     * Create a text/plain data URL from the file.  This method is used for test case generation
     * @param file
     * @return
     */
    public static String createDataURL(File file) throws IOException {
        String actual = Files.readString(Paths.get(file.getAbsolutePath()));

        return "data:," + URLEncoder.encode(actual, "utf-8");
    }

}
