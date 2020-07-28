package org.broad.igv.feature.genome;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.fasta.FastaBlockCompressedSequence;
import org.broad.igv.feature.genome.fasta.FastaIndexedSequence;
import org.broad.igv.feature.genome.fasta.FastaUtils;
import org.broad.igv.feature.gff.GFFFeatureSource;
import org.broad.igv.track.FeatureCollectionSource;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.IGV;
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
public class GenomeLoader {

    public static final long ONE_WEEK = 7 * 24 * 60 * 60 * 1000;
    private static Logger log = Logger.getLogger(GenomeLoader.class);

    static Map<String, File> localSequenceMap;

    /**
     * Define a minimal genome from a chrom.sizes file.  It is assumed (required) that the file follow the
     * UCSC naming convention  =>  [id].chrom.sizes
     *
     * @param genomePath
     * @return
     * @throws IOException
     */
    static Genome loadChromSizes(String genomePath) throws IOException {
        int firstPeriodIdx = genomePath.indexOf('.');
        String genomeId = genomePath.substring(0, firstPeriodIdx);
        List<Chromosome> chromosomes = ChromSizesParser.parse(genomePath);
        Genome newGenome = new Genome(genomeId, chromosomes);
        return newGenome;

    }

    static Genome loadGenbankFile(String genomePath) throws IOException {
        Genome newGenome;
        GenbankParser genbankParser = new GenbankParser(genomePath);
        genbankParser.readFeatures(true);

        String name = genbankParser.getLocusName();
        String chr = genbankParser.getChr();

        if (!name.equals(chr)) {
            name = name + " (" + chr + ")";
        }

        byte[] seq = genbankParser.getSequence();
        Sequence sequence = new InMemorySequence(chr, seq);
        newGenome = new Genome(chr, name, sequence, true);

        String[] aliases = genbankParser.getAliases();
        if (aliases != null) {
            List<String> aliasList = new ArrayList<String>();
            aliasList.add(chr);
            for (String a : aliases) {
                aliasList.add(a);
            }
            newGenome.addChrAliases(Arrays.asList(aliasList));
        }

        if (IGV.hasInstance() && !Globals.isHeadless()) {
            FeatureTrack geneFeatureTrack = createGeneTrack(newGenome, genbankParser.getFeatures());
            newGenome.setGeneTrack(geneFeatureTrack);
        }

        FeatureDB.addFeatures(genbankParser.getFeatures(), newGenome);

        return newGenome;
    }

    /**
     * Create a Genome from a  fasta file.
     *
     * @param genomePath
     * @return
     * @throws IOException
     */
    static Genome loadFastaFile(String genomePath) throws IOException {

        String fastaPath = null;
        String fastaIndexPath = null;
        if (genomePath.endsWith(".fai")) {
            fastaPath = genomePath.substring(0, genomePath.length() - 4);
            fastaIndexPath = genomePath;
        } else {
            fastaPath = genomePath;
            fastaIndexPath = genomePath + ".fai";
        }

        if (!FileUtils.resourceExists(fastaIndexPath)) {
            //Have to make sure we have a local copy of the fasta file
            //to index it
            if (!FileUtils.isRemote(fastaPath)) {
                fastaIndexPath = fastaPath + ".fai";
                FastaUtils.createIndexFile(fastaPath, fastaIndexPath);
            }
        }

        String id = fastaPath;
        String name;
        if (HttpUtils.isRemoteURL(fastaPath)) {
            name = Utilities.getFileNameFromURL(fastaPath);
        } else {
            File file = new File(fastaPath);
            if (!file.exists()) {
                throw new IOException(fastaPath + " does not exist, could not load genome");
            }
            name = file.getName();
        }

        FastaIndexedSequence fastaSequence = fastaPath.endsWith(".gz") ?
                new FastaBlockCompressedSequence(fastaPath) :
                new FastaIndexedSequence(fastaPath);
        Sequence sequence = new SequenceWrapper(fastaSequence);
        return new Genome(id, name, sequence, true);
    }

    static Genome loadJsonFile(String genomePath) throws IOException {

        BufferedReader reader = ParsingUtils.openBufferedReader(genomePath);
        JsonParser parser = new JsonParser();
        JsonObject json = parser.parse(reader).getAsJsonObject();

        String id = json.get("id").getAsString();
        String name = json.get("name").getAsString();
        String fastaPath = json.get("fastaURL").getAsString();
        JsonElement indexPathObject = json.get("indexURL");
        String indexPath = indexPathObject == null ? null : indexPathObject.getAsString();
        JsonElement aliasURL = json.get("aliasURL");

        fastaPath = FileUtils.getAbsolutePath(fastaPath, genomePath);
        if (indexPath != null) {
            indexPath = FileUtils.getAbsolutePath(indexPath, genomePath);
        }

        FastaIndexedSequence sequence = fastaPath.endsWith(".gz") ?
                new FastaBlockCompressedSequence(fastaPath, indexPath) :
                new FastaIndexedSequence(fastaPath, indexPath);

        ArrayList<ResourceLocator> tracks = new ArrayList<>();
        JsonArray annotations = json.getAsJsonArray("tracks");
        if (annotations == null) {
            annotations = json.getAsJsonArray("annotations");
        }
        if (annotations != null) {
            annotations.forEach((JsonElement jsonElement) -> {
                JsonObject obj = jsonElement.getAsJsonObject();
                String trackPath = obj.get("url").getAsString();
                JsonElement trackName = obj.get("name");
                JsonElement trackIndex = obj.get("indexURL");
                JsonElement indexed = obj.get("indexed");

                String trackIndexPath = null;
                if (trackPath != null) {
                    trackPath = FileUtils.getAbsolutePath(trackPath, genomePath);
                }
                if (trackIndex != null) {
                    trackIndexPath = FileUtils.getAbsolutePath(trackIndex.getAsString(), genomePath);
                }

                ResourceLocator res = new ResourceLocator(trackPath);
                if (trackName != null) res.setName(trackName.getAsString());
                if (trackIndexPath != null) res.setIndexPath(trackIndexPath);
                if (indexed != null) res.setIndexed(indexed.getAsBoolean());
                tracks.add(res);
            });
        }

        Genome newGenome = new Genome(id, name, sequence, true);
        newGenome.setAnnotationResources(tracks);
        if(aliasURL != null) {
            newGenome.addChrAliases(GenomeLoader.loadChrAliases(aliasURL.getAsString()));
        }

        return newGenome;
    }

    static Collection<Collection<String>> loadChrAliases(String path) {
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

    /**
     * Create a genome from a ".genome" file.  In addition to the reference sequence .genome files can optionally
     * specify cytobands and annotations.
     */
    static Genome loadDotGenomeFile(File archiveFile) throws IOException {

        Genome newGenome;

        if (!archiveFile.exists()) {
            throw new FileNotFoundException("Genome file: " + archiveFile.getAbsolutePath() + " does not exist.");
        }

        GenomeDescriptor genomeDescriptor = GenomeDescriptor.parseGenomeArchiveFile(archiveFile);
        final String id = genomeDescriptor.getId();
        final String displayName = genomeDescriptor.getName();
        String sequencePath = localSequenceMap.containsKey(genomeDescriptor.getId()) ?
                loadSequenceMap().get(genomeDescriptor.getId()).getAbsolutePath() :
                genomeDescriptor.getSequencePath();

        Sequence sequence;
        boolean chromosOrdered = false;
        if (sequencePath == null) {
            sequence = null;
        } else {
            if (sequencePath.endsWith(".gz")) {
                FastaBlockCompressedSequence fastaSequence = new FastaBlockCompressedSequence(sequencePath);
                sequence = new SequenceWrapper(fastaSequence);
            } else {
                FastaIndexedSequence fastaSequence = new FastaIndexedSequence(sequencePath);
                sequence = new SequenceWrapper(fastaSequence);
            }
            chromosOrdered = true;
        }

        newGenome = new Genome(id, displayName, sequence, chromosOrdered);

        if (genomeDescriptor.hasCytobands()) {
            InputStream cytobandStream = null;
            try {
                cytobandStream = genomeDescriptor.getCytoBandStream();
                BufferedReader reader = new BufferedReader(new InputStreamReader(cytobandStream));
                newGenome.setCytobands(CytoBandFileParser.loadData(reader));
            } catch (IOException ex) {
                log.warn("Error loading cytoband file", ex);
                throw new RuntimeException("Error loading cytoband file" + genomeDescriptor.cytoBandFileName);
            } finally {
                closeSilently(cytobandStream);
            }
        }


        InputStream aliasStream = null;
        try {
            aliasStream = genomeDescriptor.getChrAliasStream();
            if (aliasStream != null) {
                BufferedReader reader = new BufferedReader(new InputStreamReader(aliasStream));
                Collection<Collection<String>> aliases = loadChrAliases(reader);
                if (aliases != null) {
                    newGenome.addChrAliases(aliases);
                }
            }
        } catch (IOException e) {
            // We don't want to bomb if the alias load fails.  Just log it and proceed.
            log.error("Error loading chromosome alias table");
        } finally {
            closeSilently(aliasStream);
        }


        String geneFileName = genomeDescriptor.getGeneFileName();
        InputStream geneStream = null;
        if (geneFileName != null) {
            try {
                geneStream = genomeDescriptor.getGeneStream();
                if (geneFileName.endsWith(".gbk")) {
                    // unusual case, using genbank file for features within a .genome file
                    GenbankParser genbankParser = new GenbankParser();
                    genbankParser.readFeatures(geneStream, false);
                    FeatureTrack geneFeatureTrack = createGeneTrack(newGenome, genbankParser.getFeatures());
                    newGenome.setGeneTrack(geneFeatureTrack);
                } else {
                    BufferedReader reader = new BufferedReader(new InputStreamReader(geneStream));
                    FeatureTrack geneFeatureTrack = createGeneTrack(newGenome, reader,
                            geneFileName, genomeDescriptor.getGeneTrackName(),
                            genomeDescriptor.getUrl());

                    newGenome.setGeneTrack(geneFeatureTrack);
                }
            } finally {
                closeSilently(geneStream);
            }
        }

        genomeDescriptor.close();

        return newGenome;
    }

    private static Collection<Collection<String>> loadChrAliases(BufferedReader br) throws IOException {
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
     * @param reader        a reader for the gene (annotation) file.
     * @param genome
     * @param geneFileName
     * @param geneTrackName
     */
    public static FeatureTrack createGeneTrack(Genome genome, BufferedReader reader, String geneFileName, String geneTrackName,
                                               String annotationURL) {

        FeatureDB.clearFeatures();
        FeatureTrack geneFeatureTrack = null;

        if (reader != null) {
            FeatureParser parser;
            if (geneFileName.endsWith(".embl")) {
                parser = new EmblFeatureTableParser();
            } else if (GFFFeatureSource.isGFF(geneFileName)) {
                parser = new GFFParser();
            } else {
                parser = AbstractFeatureParser.getInstanceFor(new ResourceLocator(geneFileName), genome);
            }
            if (parser == null) {
                MessageUtils.showMessage("ERROR: Unrecognized annotation file format: " + geneFileName +
                        "<br>Annotations for genome: " + genome.getId() + " will not be loaded.");
            } else {
                List<htsjdk.tribble.Feature> genes = parser.loadFeatures(reader, genome);
                String name = geneTrackName;
                if (name == null) name = "Genes";

                String id = genome.getId() + "_genes";
                geneFeatureTrack = new FeatureTrack(id, name, new FeatureCollectionSource(genes, genome));
                geneFeatureTrack.setMinimumHeight(5);
                geneFeatureTrack.setHeight(35);
                geneFeatureTrack.setTrackType(TrackType.GENE);
                geneFeatureTrack.setColor(Color.BLUE.darker());
                TrackProperties props = parser.getTrackProperties();
                if (props != null) {
                    geneFeatureTrack.setProperties(parser.getTrackProperties());
                }
                geneFeatureTrack.setUrl(annotationURL);
            }
        }
        return geneFeatureTrack;
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

    static Map<String, File> loadSequenceMap() {

        File sequenceFile = new File(DirectoryManager.getGenomeCacheDirectory(), GenomeManager.SEQUENCE_MAP_FILE);
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

    private static void closeSilently(InputStream stream) {
        if(stream != null) {
            try {
                stream.close();
            } catch (IOException e) {
                log.error("Error closing stream", e);
            }
        }
    }

    private static void closeSilently(Reader reader) {
        if(reader != null) {
            try {
                reader.close();
            } catch (IOException e) {
                log.error("Error closing reader", e);
            }
        }
    }
}
