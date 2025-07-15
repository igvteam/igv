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
import org.broad.igv.util.FileUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Field;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.*;

/**
 * Class of static functions for managing genome downloads
 */
public class GenomeDownloadUtils {

    private static Logger log = LogManager.getLogger(GenomeDownloadUtils.class);

    public static boolean isAnnotationsDownloadable(String path) {
        return path.endsWith(".json");
    }

    public static boolean isSequenceDownloadable(String genomePath) {
        if (genomePath != null && genomePath.endsWith(".json")) {
            try {
                String jsonString = FileUtils.getContents(genomePath);
                GenomeConfig genomeConfig = GenomeConfig.fromJson(jsonString);
                String sequenceURL = genomeConfig.twoBitURL;
                if (sequenceURL == null) {
                    sequenceURL = genomeConfig.fastaURL;
                }
                return isRemoteURL(sequenceURL); // && !disallowedBuckets.stream().anyMatch(sequenceURL::contains);

            } catch (IOException e) {
                log.error("Error fetching genome json " + genomePath);
            }
        }
        return false;

    }

    public static File downloadGenome(GenomeConfig genomeConfig, boolean downloadSequence, boolean downloadAnnotations) throws IOException {

        GenomeConfig config = genomeConfig.copy();

        // Create a directory for the data files (sequence and annotations)
        final File genomeDirectory = DirectoryManager.getGenomeCacheDirectory();
        File dataDirectory = null;

        if(downloadSequence || downloadAnnotations) {
            dataDirectory = new File(genomeDirectory, config.id);
            if (dataDirectory.exists()) {
                if (!dataDirectory.isDirectory()) {
                    throw new RuntimeException("Error downloading genome. " + dataDirectory.getAbsolutePath() + " exists and is not a directory.");
                }
            } else {
                if (!dataDirectory.mkdir()) {
                    throw new RuntimeException("Error downloading genome.  Could not create directory: " + dataDirectory.getAbsolutePath());
                }
            }
        }

        String relativeDataDirectory = config.id + "/";

        if (config.twoBitURL != null) {

            URL url = new URL(config.twoBitURL);
            File localFile;
            if (downloadSequence) {
                localFile = download(url, dataDirectory);
            } else {
                localFile = constructLocalFile(url, dataDirectory);  // It might be there from previous downloads
            }
            if (localFile.exists()) {
                config.twoBitURL =(relativeDataDirectory + localFile.getName());
                // Null out urls not needed for .2bit sequences.
                config.fastaURL = (null);
                config.indexURL = (null);
                config.gziIndexURL = (null);
                config.twoBitBptURL =(null);  // Not needed for local .2bit files
            }
        } else if (config.fastaURL != null) {

            String[] fastaFields = {"fastaURL", "indexURL", "gziIndexURL", "compressedIndexURL"};
            downloadAndUpdateConfig(downloadSequence, fastaFields, config, dataDirectory, relativeDataDirectory, true);

        } else {
            throw new RuntimeException("Sequence for genome " + config.getName() + " is not downloadable.");
        }

        String[] annotationFields = {"chromAliasBbURL", "cytobandURL", "cytobandBbURL", "aliasURL", "chromSizesURL"};
        downloadAndUpdateConfig(downloadAnnotations, annotationFields, config, dataDirectory, relativeDataDirectory, false);

        List<TrackConfig> trackConfigs = config.getTrackConfigs();
        if (trackConfigs != null) {
            String[] trackFields = {"url", "indexURL", "trixURL"};
            for (TrackConfig trackConfig : trackConfigs) {
                downloadAndUpdateConfig(downloadAnnotations, trackFields, trackConfig, dataDirectory, relativeDataDirectory, false);
            }

        }

        return saveLocalGenome(config);
    }

    private static void downloadAndUpdateConfig(
            boolean downloadData,
            String[] fields,
            Object config,
            File dataDirectory,
            String relativeDataDirectory,
            boolean useLocalCache) {
        URL url;
        File localFile = null;
        for (String f : fields) {
            try {
                Field field = config.getClass().getDeclaredField(f);
                field.setAccessible(true);
                Object urlField = field.get(config);
                if (urlField != null) {
                    url = new URL(urlField.toString());
                    if (downloadData) {
                        localFile = download(url, dataDirectory);
                    } else if (useLocalCache) {
                        localFile = constructLocalFile(url, dataDirectory);   // It might be there from previous downloads
                    }
                    if (localFile != null && localFile.exists()) {
                        field.set(config, relativeDataDirectory + localFile.getName());
                    }
                }
            } catch (Exception e) {
                log.error(e);
            }
        }
    }

    /**
     * Save the genome definition as a json file in the genome directory.
     *
     * @param genomeConfig
     * @return
     * @throws IOException
     */
    public static File saveLocalGenome(GenomeConfig genomeConfig) throws IOException {
        String id = genomeConfig.id;
        String sanitizedId = id.replaceAll("[^a-zA-Z0-9-_]", "_");
        File localFile = new File(DirectoryManager.getGenomeCacheDirectory(), sanitizedId + ".json");
        log.info("Saving " + localFile.getAbsolutePath());
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(localFile))) {
            Gson gson = new GsonBuilder().setPrettyPrinting().create();
            gson.toJson(genomeConfig, writer);
        }
        return localFile;
    }

    public static File constructLocalFile(URL url, File directory) {
        String path = url.getPath();
        int idx = path.lastIndexOf('/');
        String filename = idx < 0 ? path : path.substring(idx + 1);
        return new File(directory, filename);
    }

    private static File download(URL url, File directory) throws MalformedURLException {

        File localFile = constructLocalFile(url, directory);

        boolean success = Downloader.download(url, localFile, IGV.getInstance().getMainFrame());
        if (!success) {
            throw new RuntimeException("Download canceled");
        }
        return localFile;
    }

    private static boolean isRemoteURL(String url) {
        return url.startsWith("http://") || url.startsWith("https://");
    }

    static Set<String> disallowedBuckets = new HashSet<>(Arrays.asList(
            "igv.org.genomes",
            "igv-genepattern-org",
            "igv.broadinstitute.org"));

}
