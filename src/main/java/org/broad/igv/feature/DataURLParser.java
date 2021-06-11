package org.broad.igv.feature;

import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.track.FeatureSource;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;

import java.io.*;
import java.net.URLDecoder;
import java.net.URLEncoder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Iterator;
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

        String fn = "ignore." + format.toLowerCase();
        ResourceLocator locator = new ResourceLocator(fn);   // "fake" location
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
