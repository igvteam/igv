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
        File dataDirectory = new File(genomeDirectory, config.getId());
        if (dataDirectory.exists()) {
            if (!dataDirectory.isDirectory()) {
                throw new RuntimeException("Error downloading genome. " + dataDirectory.getAbsolutePath() + " exists and is not a directory.");
            }
        } else {
            if (!dataDirectory.mkdir()) {
                throw new RuntimeException("Error downloading genome.  Could not create directory: " + dataDirectory.getAbsolutePath());
            }
        }

        String relativeDataDirectory = config.getId() + "/";
        File localFile;
        if (downloadSequence) {

            if (config.getTwoBitURL() != null) {

                localFile = download(new URL(config.getTwoBitURL()), dataDirectory);
                config.setTwoBitURL(relativeDataDirectory + localFile.getName());

                // Null out urls not needed for .2bit seequences.
                config.setFastaURL(null);
                config.setIndexURL(null);
                config.setGziIndexURL(null);

                // The optional btree index.  This is probably not needed for local file configurations.
                if (config.getTwoBitBptURL() != null) {
                    localFile = download(new URL(config.getTwoBitBptURL()), dataDirectory);
                    config.setTwoBitBptURL(relativeDataDirectory + localFile.getName());
                }
            } else if (config.getFastaURL() != null) {

                localFile = download(new URL(config.getFastaURL()), dataDirectory);
                config.setFastaURL(relativeDataDirectory + localFile.getName());
                if (config.getIndexURL() != null) {
                    localFile = download(new URL(config.getIndexURL()), dataDirectory);
                    config.setIndexURL(relativeDataDirectory + localFile.getName());
                }
                if (config.getGziIndexURL() != null) {
                    localFile = download(new URL(config.getGziIndexURL()), dataDirectory);
                    config.setGziIndexURL(relativeDataDirectory + localFile.getName());
                }

            } else {
                throw new RuntimeException("Sequence for genome " + config.getName() + " is not downloadable.");
            }
        }
        if (downloadAnnotations) {

            if (config.getChromSizesURL() != null) {
                localFile = download(new URL(config.getChromSizesURL()), dataDirectory);
                config.setChromSizesURL(relativeDataDirectory + localFile.getName());
            }

            if (config.getCytobandBbURL() != null) {
                localFile = download(new URL(config.getCytobandBbURL()), dataDirectory);
                config.setCytobandBbURL(relativeDataDirectory + localFile.getName());

            } else if (config.getCytobandURL() != null) {
                localFile = download(new URL(config.getCytobandURL()), dataDirectory);
                config.setCytobandURL(relativeDataDirectory + localFile.getName());
            }

            if (config.getChromAliasBbURL() != null) {
                localFile = download(new URL(config.getChromAliasBbURL()), dataDirectory);
                config.setChromAliasBbURL(relativeDataDirectory + localFile.getName());

            } else if (config.getAliasURL() != null) {
                localFile = download(new URL(config.getAliasURL()), dataDirectory);
                config.setAliasURL(relativeDataDirectory + localFile.getName());
            }

            List<TrackConfig> trackConfigs = config.getTrackConfigs();

            if (trackConfigs != null) {
                for (TrackConfig trackConfig : trackConfigs) {
                    if (trackConfig.getUrl() != null) {
                        localFile = download(new URL(trackConfig.getUrl()), dataDirectory);
                        trackConfig.setUrl(relativeDataDirectory + localFile.getName());
                    }
                    if (trackConfig.getIndexURL() != null) {
                        localFile = download(new URL(trackConfig.getIndexURL()), dataDirectory);
                        trackConfig.setIndexURL(relativeDataDirectory + localFile.getName());
                    }
                    if(trackConfig.getTrixURL() != null) {
                        localFile = download(new URL(trackConfig.getTrixURL()), dataDirectory);
                        trackConfig.setTrixURL(relativeDataDirectory + localFile.getName());
                    }
                }

            }
        }
        File localGenomeFile = new File(genomeDirectory, config.getId() + ".json");
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
