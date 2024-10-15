package org.broad.igv.feature.genome;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import org.broad.igv.DirectoryManager;
import org.broad.igv.feature.genome.load.GenomeConfig;
import org.broad.igv.feature.genome.load.TrackConfig;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.download.Downloader;
import org.broad.igv.util.HttpUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.List;

/**
 * Class of static functions for managing genome downloads
 */
public class GenomeDownloadUtils {

    private static Logger log = LogManager.getLogger(GenomeDownloadUtils.class);

    public static boolean isAnnotationsDownloadable(GenomeListItem item) {
        return item.getPath().endsWith(".json");
    }

    public static boolean isSequenceDownloadable(GenomeListItem item) {
        if (item.getPath().endsWith(".json")) {
            try {
                String jsonString = HttpUtils.getInstance().getContentsAsJSON(new URL(item.getPath()));
                return jsonString.contains("twoBitURL");
            } catch (IOException e) {
                log.error("Error fetching genome json " + item.getPath());
            }
        }
        return false;
    }

    public static File downloadGenome(GenomeConfig config, boolean downloadSequence, boolean downloadAnnotations) throws IOException {

        // Create a directory for the data files (sequence and annotations)
        final File genomeDirectory = DirectoryManager.getGenomeCacheDirectory();
        File dataDirectory = new File(genomeDirectory, config.id);
        if (dataDirectory.exists()) {
            if (!dataDirectory.isDirectory()) {
                throw new RuntimeException("Error downloading genome. " + dataDirectory.getAbsolutePath() + " exists and is not a directory.");
            }
        } else {
            if (!dataDirectory.mkdir()) {
                throw new RuntimeException("Error downloading genome.  Could not create directory: " + dataDirectory.getAbsolutePath());
            }
        }

        String relativeDataDirectory = config.id + "/";
        File localFile;
        if (downloadSequence) {

            if (config.twoBitURL != null) {

                localFile = download(new URL(config.twoBitURL), dataDirectory);
                config.twoBitURL = relativeDataDirectory + localFile.getName();

                // Null out urls not needed for .2bit seequences.
                config.fastaURL = null;
                config.indexURL = null;
                config.gziIndexURL = null;

                // The optional btree index.  This is probably not needed for local file configurations.
                if (config.twoBitBptURL != null) {
                    localFile = download(new URL(config.twoBitBptURL), dataDirectory);
                    config.twoBitBptURL = relativeDataDirectory + localFile.getName();
                }
            } else if (config.fastaURL != null) {

                localFile = download(new URL(config.fastaURL), dataDirectory);
                config.fastaURL = relativeDataDirectory + localFile.getName();
                if (config.indexURL != null) {
                    localFile = download(new URL(config.indexURL), dataDirectory);
                    config.indexURL = relativeDataDirectory + localFile.getName();
                }
                if (config.gziIndexURL != null) {
                    localFile = download(new URL(config.gziIndexURL), dataDirectory);
                    config.gziIndexURL = relativeDataDirectory + localFile.getName();
                }

            } else {
                throw new RuntimeException("Sequence for genome " + config.name + " is not downloadable.");
            }
        }
        if (downloadAnnotations) {

            if (config.chromSizesURL != null) {
                localFile = download(new URL(config.chromSizesURL), dataDirectory);
                config.chromSizesURL = relativeDataDirectory + localFile.getName();
            }

            if (config.cytobandBbURL != null) {
                localFile = download(new URL(config.cytobandBbURL), dataDirectory);
                config.cytobandBbURL = relativeDataDirectory + localFile.getName();

            } else if (config.cytobandURL != null) {
                localFile = download(new URL(config.cytobandURL), dataDirectory);
                config.cytobandURL = relativeDataDirectory + localFile.getName();
            }

            if (config.chromAliasBbURL != null) {
                localFile = download(new URL(config.chromAliasBbURL), dataDirectory);
                config.chromAliasBbURL = relativeDataDirectory + localFile.getName();

            } else if (config.aliasURL != null) {
                localFile = download(new URL(config.aliasURL), dataDirectory);
                config.aliasURL = relativeDataDirectory + localFile.getName();
            }

            List<TrackConfig> trackConfigs = config.tracks;
            if (trackConfigs == null) trackConfigs = config.annotations;   // alias

            if (trackConfigs != null) {
                for (TrackConfig trackConfig : trackConfigs) {
                    if (trackConfig.url != null) {
                        localFile = download(new URL(trackConfig.url), dataDirectory);
                        trackConfig.url = relativeDataDirectory + localFile.getName();
                    }
                    if (trackConfig.indexURL != null) {
                        localFile = download(new URL(trackConfig.indexURL), dataDirectory);
                        trackConfig.indexURL = relativeDataDirectory + localFile.getName();
                    }
                }

            }
        }
        File localGenomeFile = new File(genomeDirectory, config.id + ".json");
        saveLocalGenome(config, localGenomeFile);
        return localGenomeFile;

    }

    private static void saveLocalGenome(GenomeConfig genomeConfig, File localFile) throws IOException {
        log.info("Saving " + localFile.getAbsolutePath());
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(localFile))) {
            Gson gson = new GsonBuilder().setPrettyPrinting().create();
            gson.toJson(genomeConfig, writer);
        }
    }

    private static File download(URL url, File directory) throws MalformedURLException {
        String path = url.getPath();
        int idx = path.lastIndexOf('/');
        String filename = idx < 0 ? path : path.substring(idx + 1);
        File localFile = new File(directory, filename);
        boolean success = Downloader.download(url, localFile, IGV.getInstance().getMainFrame());
        if (!success) {
            throw new RuntimeException("Download canceled");
        }
        return localFile;
    }

}
