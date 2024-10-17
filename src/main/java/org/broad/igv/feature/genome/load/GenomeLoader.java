package org.broad.igv.feature.genome.load;

import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.genome.DotGenomeUtils;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeException;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.track.FeatureCollectionSource;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.Utilities;

import java.awt.*;
import java.io.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Collection of static load methods for various genome definition formats
 * - genome files (legacy)
 * - genbank
 * - json
 * - fasta
 * - chrom sizes (minimal, no reference sequence)
 */
abstract public class GenomeLoader {

    private static Logger log = LogManager.getLogger(GenomeLoader.class);

    public static Map<String, File> localSequenceMap;

    public static GenomeLoader getLoader(String genomePath) throws IOException {
        if (genomePath.endsWith(".genome")) {
            File archiveFile = DotGenomeUtils.getDotGenomeFile(genomePath);
            if (!archiveFile.exists()) {
                MessageUtils.showMessage("File not found: " + archiveFile.getAbsolutePath());
                return null;    // Happens if genome download was canceled.
            }
            return new DotGenomeLoader(archiveFile);
        } else if (genomePath.endsWith(".gbk") || genomePath.endsWith(".gb")) {
            return new GenbankLoader(genomePath);
        } else if (genomePath.endsWith(".chrom.sizes")) {
            return new ChromsizesLoader(genomePath);
        } else if (genomePath.endsWith(".json")) {
            return new JsonGenomeLoader(genomePath);
        } else if (genomePath.endsWith(".2bit")) {
            GenomeConfig config = new GenomeConfig();
            config.setTwoBitURL(genomePath);
            config.setId(genomePath);
            config.setName((HttpUtils.isRemoteURL(genomePath)) ?
                    Utilities.getFileNameFromURL(genomePath) :
                    (new File(genomePath)).getName());
            return new GenomeObjectLoader(config);
        } else if (genomePath.endsWith("hub.txt")) {
            return new HubGenomeLoader(genomePath);
        } else {
            // Assume a fasta or 2bit file file
            if (genomePath.endsWith(Globals.GZIP_FILE_EXTENSION)) {
                String gziPath = genomePath + ".gzi";
                String faiPath = genomePath + ".fai";
                if (!(FileUtils.resourceExists(gziPath) && FileUtils.resourceExists(faiPath))) {
                    throw new GenomeException("IGV cannot readed gzipped fasta files.");
                }
            }
            return new FastaGenomeLoader(genomePath);
        }
    }

    abstract public Genome loadGenome() throws IOException;


    /**
     * Create an annotation track for the genome from a supplied list of features
     *
     * @param genome
     * @param features
     */
    public static FeatureTrack createGeneTrack(Genome genome, List<htsjdk.tribble.Feature> features) {

        FeatureDB.clearFeatures();
        FeatureTrack geneFeatureTrack = null;
        String name = "Annotations";

        String id = genome.getId() + "_genes";
        geneFeatureTrack = new FeatureTrack(id, name, new FeatureCollectionSource(features, genome));
        geneFeatureTrack.setMinimumHeight(5);
        geneFeatureTrack.setHeight(35);
        //geneFeatureTrack.setRendererClass(GeneRenderer.class);
        geneFeatureTrack.setColor(Color.BLUE.darker());

        return geneFeatureTrack;
    }

    public static Map<String, File> loadSequenceMap() {

        File sequenceFile = new File(DirectoryManager.getGenomeCacheDirectory(), GenomeDescriptor.SEQUENCE_MAP_FILE);
        localSequenceMap = new HashMap<>();
        if (sequenceFile.exists()) {
            BufferedReader br = null;
            try {
                br = new BufferedReader(new FileReader(sequenceFile));
                String nextLine;
                while ((nextLine = br.readLine()) != null) {
                    String[] tokens = nextLine.split("\t");
                    if (tokens.length > 1) {
                        File file = new File(tokens[1]);
                        if (file.exists()) {
                            localSequenceMap.put(tokens[0], file);
                        } else {
                            log.warn("Sequence file not found: " + file.getAbsolutePath());
                        }
                    }
                }
            } catch (IOException e) {
                log.error("Error loading sequence map file", e);
            } finally {
                closeSilently(br);
            }
        }
        return localSequenceMap;
    }

    static void closeSilently(InputStream stream) {
        if (stream != null) {
            try {
                stream.close();
            } catch (IOException e) {
                log.error("Error closing stream", e);
            }
        }
    }

    static void closeSilently(Reader reader) {
        if (reader != null) {
            try {
                reader.close();
            } catch (IOException e) {
                log.error("Error closing reader", e);
            }
        }
    }
}
