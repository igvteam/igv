package org.igv.feature.genome.load;

import org.igv.feature.genome.Genome;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.util.FileUtils;
import org.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.util.List;

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
        return genome;
    }

    public GenomeConfig loadGenomeConfig() throws IOException {

        try (InputStream is = ParsingUtils.openInputStream(genomePath)) {
            String jsonString = ParsingUtils.readContentsFromStream(is);
            GenomeConfig genomeConfig = GenomeConfig.fromJson(jsonString);
            fixPaths(genomeConfig);
            return genomeConfig;
        }
    }

    /**
     * JSON paths/urls can be relative to the path to the genome json file.  This method converts them to absolulte
     * paths/urls.
     *
     * @param config
     * @return
     */
    private GenomeConfig fixPaths(GenomeConfig config) {

        if (config.chromAliasBbURL != null) {
            config.chromAliasBbURL = FileUtils.getAbsolutePath(config.chromAliasBbURL, genomePath);
        }
        if (config.twoBitURL != null) {
            config.twoBitURL = FileUtils.getAbsolutePath(config.twoBitURL, genomePath);
        }
        if (config.twoBitBptURL != null) {
            config.twoBitBptURL = FileUtils.getAbsolutePath(config.twoBitBptURL, genomePath);
        }
        if (config.cytobandBbURL != null) {
            config.cytobandBbURL = FileUtils.getAbsolutePath(config.cytobandBbURL, genomePath);
        }
        if (config.chromSizesURL != null) {
            config.chromSizesURL = FileUtils.getAbsolutePath(config.chromSizesURL, genomePath);
        }
        if (config.fastaURL != null) {
            config.fastaURL = FileUtils.getAbsolutePath(config.fastaURL, genomePath);
        }
        if (config.indexURL != null) {
            config.indexURL = FileUtils.getAbsolutePath(config.indexURL, genomePath);
        }
        if (config.gziIndexURL != null) {
            config.gziIndexURL = FileUtils.getAbsolutePath(config.gziIndexURL, genomePath);
        }
        if (config.cytobandURL != null) {
            config.cytobandURL = FileUtils.getAbsolutePath(config.cytobandURL, genomePath);
        }
        if (config.aliasURL != null) {
            config.aliasURL = FileUtils.getAbsolutePath(config.aliasURL, genomePath);
        }
        if (config.chromAliasBbURL != null) {
            config.chromAliasBbURL = FileUtils.getAbsolutePath(config.chromAliasBbURL, genomePath);
        }
        if (config.infoURL != null) {
            config.infoURL = FileUtils.getAbsolutePath(config.infoURL, genomePath);
        }

        List<TrackConfig> trackConfigs = config.getTrackConfigs();
        if (trackConfigs != null) {
            trackConfigs.forEach((TrackConfig trackConfig) -> {
                if (trackConfig.url != null) {
                    trackConfig.url = (FileUtils.getAbsolutePath(trackConfig.url, genomePath));
                }
                if (trackConfig.indexURL != null) {
                    trackConfig.indexURL = (FileUtils.getAbsolutePath(trackConfig.indexURL, genomePath));
                }
            });
        }
        return config;
    }


    public GenomeDescriptor loadDescriptor() throws IOException {
        BufferedReader reader = null;
        try {
            reader = ParsingUtils.openBufferedReader(genomePath);
            StringBuilder sb = new StringBuilder();
            String line;
            while ((line = reader.readLine()) != null) {
                sb.append(line);
            }
            org.json.JSONObject json = new org.json.JSONObject(sb.toString());
            String id = json.optString("id");
            String name = json.optString("name");
            String fastaPath = json.optString("fastaURL");
            return new GenomeDescriptor(id, name, fastaPath);
        } finally {
            if (reader != null) {
                reader.close();
            }
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

    }

}
