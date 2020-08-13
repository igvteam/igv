package org.broad.igv.feature.genome.load;

import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.*;
import org.broad.igv.track.FeatureCollectionSource;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.*;

import java.awt.*;
import java.io.*;
import java.util.List;
import java.util.*;

/**
 * Collection of static load methods for various genome definition formats
 * - genome files (legacy)
 * - genbank
 * - json
 * - fasta
 * - chrom sizes (minimal, no reference sequence)
 */
abstract public class GenomeLoader {

    public static final long ONE_WEEK = 7 * 24 * 60 * 60 * 1000;
    private static Logger log = Logger.getLogger(GenomeLoader.class);

    public static Map<String, File> localSequenceMap;

    public static GenomeLoader getLoader(String genomePath) throws IOException {
        if (genomePath.endsWith(".genome")) {
            File archiveFile = GenomeManager.getGenomeFile(genomePath);
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
        } else {
            // Assume a fasta file
            if (genomePath.endsWith(Globals.GZIP_FILE_EXTENSION)) {
                String gziPath = genomePath + ".gzi";
                String faiPath = genomePath + ".fai";
                if (!(FileUtils.resourceExists(gziPath) && FileUtils.resourceExists(faiPath))) {
                    throw new GenomeException("IGV cannot readed gzipped fasta files.");
                }
            }
            if (!FileUtils.isRemote(genomePath)) {
                if (!(new File(genomePath)).exists()) {
                    throw new GenomeException("Cannot locate genome: " + genomePath);
                }
            }
            return new FastaGenomeLoader(genomePath);
        }
    }

    abstract public Genome loadGenome() throws IOException;

    public static Collection<Collection<String>> loadChrAliases(String path) {
        File aliasFile = new File(path);
        if (aliasFile.exists()) {
            BufferedReader br = null;
            try {
                br = new BufferedReader(new FileReader(aliasFile));
                return loadChrAliases(br);
            } catch (IOException e) {
                log.error("Error loading chr alias table", e);
                MessageUtils.showMessage("<html>Error loading chromosome alias table.  Aliases will not be avaliable<br>" +
                        e.toString());
            } finally {
                closeSilently(br);
            }
        }
        return null;
    }

    static Collection<Collection<String>> loadChrAliases(BufferedReader br) throws IOException {
        String nextLine = "";
        Collection<Collection<String>> synonymList = new ArrayList<Collection<String>>();
        while ((nextLine = br.readLine()) != null) {
            String[] tokens = nextLine.split("\t");
            if (tokens.length > 1) {
                Collection<String> synonyms = new ArrayList<String>();
                for (String t : tokens) {
                    String syn = t.trim();
                    if (t.length() > 0) synonyms.add(syn.trim());
                }
                synonymList.add(synonyms);
            }
        }
        return synonymList;
    }

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
                            log.info("Sequence file not found: " + file.getAbsolutePath());
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
