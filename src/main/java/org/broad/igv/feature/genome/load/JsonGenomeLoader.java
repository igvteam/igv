package org.broad.igv.feature.genome.load;

import com.google.gson.Gson;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureReader;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.IGVNamedFeature;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.track.TribbleFeatureSource;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

public class JsonGenomeLoader extends GenomeLoader {

    private static Logger log = LogManager.getLogger(JsonGenomeLoader.class);

    private String genomePath;

    public JsonGenomeLoader(String genomePath) {
        this.genomePath = genomePath;
    }

    @Override
    public Genome loadGenome() throws IOException {


        GenomeConfig genomeConfig = loadGenomeConfig();

        Genome genome = new Genome(genomeConfig);

        // Load liftover "chain" files.  This enables navigating by coordinates of another genome.
        // Not a common option.
//            JsonElement chains = config.chains;
//            if (chains != null) {
//                Map<String, Liftover> liftoverMap = new HashMap<>();
//                JsonObject chainsObj = chains.getAsJsonObject();
//                for (Map.Entry<String, JsonElement> entry : chainsObj.entrySet()) {
//                    String chainsPath = FileUtils.getAbsolutePath(entry.getValue().getAsString(), genomePath);
//                    liftoverMap.put(entry.getKey(), Liftover.load(chainsPath));
//
//                }
//                newGenome.setLiftoverMap(liftoverMap);
//            }

        return genome;


    }

    public GenomeConfig loadGenomeConfig() throws IOException {

        try (InputStream is = ParsingUtils.openInputStream(genomePath)) {

            String jsonString = ParsingUtils.readContentsFromStream(is);

            if (jsonString.contains("chromosomeOrder")) {
                jsonString = fixChromosomeOrder(jsonString);
            }

            GenomeConfig genomeConfig = GenomeConfig.fromJson(jsonString);

            fixPaths(genomeConfig);

            return genomeConfig;

        }
    }

    /**
     * Fix deprecated form of chromosome order (comma delimited list of strings)
     *
     * @param jsonString
     * @return
     */
    private String fixChromosomeOrder(String jsonString) {
        Map obj = (new Gson()).fromJson(jsonString, Map.class);
        Object chromosomeOrder = obj.get("chromosomeOrder");
        if (chromosomeOrder != null) {
            if (chromosomeOrder instanceof String) {
                obj.put("chromosomeOrder", Arrays.stream(((String) chromosomeOrder).split(",")).map(c -> c.trim()).toArray());
            }
        }
        return (new Gson()).toJson(obj);

    }

    /**
     * JSON paths/urls can be relative to the path to the genome json file.  This method converts them to absolulte
     * paths/urls.
     *
     * @param config
     * @return
     */
    private GenomeConfig fixPaths(GenomeConfig config) {

        if (config.getChromAliasBbURL() != null) {
            config.setChromAliasBbURL(FileUtils.getAbsolutePath(config.getChromAliasBbURL(), genomePath));
        }
        if (config.getTwoBitURL() != null) {
            config.setTwoBitURL(FileUtils.getAbsolutePath(config.getTwoBitURL(), genomePath));
        }
        if (config.getTwoBitBptURL() != null) {
            config.setTwoBitBptURL(FileUtils.getAbsolutePath(config.getTwoBitBptURL(), genomePath));
        }
        if (config.getCytobandBbURL() != null) {
            config.setCytobandBbURL(FileUtils.getAbsolutePath(config.getCytobandBbURL(), genomePath));
        }
        if (config.getChromSizesURL() != null) {
            config.setChromSizesURL(FileUtils.getAbsolutePath(config.getChromSizesURL(), genomePath));
        }
        if (config.getFastaURL() != null) {
            config.setFastaURL(FileUtils.getAbsolutePath(config.getFastaURL(), genomePath));
        }
        if (config.getIndexURL() != null) {
            config.setIndexURL(FileUtils.getAbsolutePath(config.getIndexURL(), genomePath));
        }
        if (config.getGziIndexURL() != null) {
            config.setGziIndexURL(FileUtils.getAbsolutePath(config.getGziIndexURL(), genomePath));
        }
        if (config.getCytobandURL() != null) {
            config.setCytobandURL(FileUtils.getAbsolutePath(config.getCytobandURL(), genomePath));
        }
        if (config.getAliasURL() != null) {
            config.setAliasURL(FileUtils.getAbsolutePath(config.getAliasURL(), genomePath));
        }
        if (config.getChromAliasBbURL() != null) {
            config.setChromAliasBbURL(FileUtils.getAbsolutePath(config.getChromAliasBbURL(), genomePath));
        }
        List<TrackConfig> trackConfigs = config.getTrackConfigs();
        if (trackConfigs != null) {
            trackConfigs.forEach((TrackConfig trackConfig) -> {
                if (trackConfig.getUrl() != null) {
                    trackConfig.setUrl(FileUtils.getAbsolutePath(trackConfig.getUrl(), genomePath));
                }
                if (trackConfig.getIndexURL() != null) {
                    trackConfig.setIndexURL(FileUtils.getAbsolutePath(trackConfig.getIndexURL(), genomePath));
                }
            });
        }
        return config;
    }


    public GenomeDescriptor loadDescriptor() throws IOException {
        BufferedReader reader = null;
        try {
            reader = ParsingUtils.openBufferedReader(genomePath);
            JsonParser parser = new JsonParser();
            JsonObject json = parser.parse(reader).getAsJsonObject();
            String id = json.get("id").getAsString();
            String name = json.get("name").getAsString();
            String fastaPath = json.get("fastaURL").getAsString();
            return new GenomeDescriptor(id, name, fastaPath);
        } finally {
            reader.close();
        }
    }

    private void addToFeatureDB(List<ResourceLocator> locators, Genome genome) {
        for (ResourceLocator locator : locators) {
            try {
                FeatureReader featureReader = TribbleFeatureSource.getBasicReader(locator, genome);
                CloseableTribbleIterator<Feature> iter = featureReader.iterator();
                while (iter.hasNext()) {
                    Feature f = iter.next();
                    if (f instanceof IGVNamedFeature) {
                        FeatureDB.addFeature((IGVNamedFeature) f, genome);
                    }
                }
            } catch (IOException e) {
                log.error("Error loading " + locator.getPath());
            }
        }
    }

    private String stripQuotes(String str) {
        if (str.startsWith("\"")) {
            return str.substring(1, str.length() - 1);  // Assume also ends with
        } else {
            return str;
        }
    }

    public static class GenomeDescriptor {
        String id;
        String name;
        String fastaURL;

        public GenomeDescriptor(String id, String name, String fastaURL) {
            this.id = id;
            this.name = name;
            this.fastaURL = fastaURL;
        }

        public String getId() {
            return id;
        }

        public String getName() {
            return name;
        }

        public String getFastaURL() {
            return fastaURL;
        }
    }

}
