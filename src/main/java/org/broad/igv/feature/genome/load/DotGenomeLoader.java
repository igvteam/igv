package org.broad.igv.feature.genome.load;

import htsjdk.tribble.Feature;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.gff.GFFFeatureSource;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.track.FeatureCollectionSource;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.io.*;
import java.util.List;


/**
 * Load for the legacy ".genome" format.   This format id deprecated.
 */
public class DotGenomeLoader extends GenomeLoader {

    private static Logger log = LogManager.getLogger(DotGenomeLoader.class);

    File archiveFile;

    public DotGenomeLoader(File archiveFile) {
        this.archiveFile = archiveFile;
    }

    /**
     * @param reader        a reader for the gene (annotation) file.
     * @param genome
     * @param geneFileName
     * @param geneTrackName
     */
    private static FeatureTrack createGeneTrack(Genome genome, BufferedReader reader, String geneFileName, String geneTrackName,
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
                List<Feature> genes = parser.loadFeatures(reader, genome);
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
                geneFeatureTrack.setFeatureInfoURL(annotationURL);
            }
        }
        return geneFeatureTrack;
    }

    /**
     * Create a genome from a ".genome" file.  In addition to the reference sequence .genome files can optionally
     * specify cytobands and annotations.
     */

    @Override
    public Genome loadGenome() throws IOException {

        GenomeConfig config = new GenomeConfig();

        GenomeDescriptor genomeDescriptor = GenomeDescriptor.parseGenomeArchiveFile(archiveFile);
        final String id = genomeDescriptor.getId();
        final String displayName = genomeDescriptor.getName();
        String sequencePath = localSequenceMap.containsKey(genomeDescriptor.getId()) ?
                loadSequenceMap().get(genomeDescriptor.getId()).getAbsolutePath() :
                genomeDescriptor.getSequencePath();

        if (sequencePath == null) {
            // TODO -- is this an error?
        } else {
            config.setFastaURL(sequencePath);
            config.setIndexURL(sequencePath + ".fai");
            if (sequencePath.endsWith(".gz")) {
                config.setGziIndexURL(sequencePath + ".gzi");
            }
        }

        config.setId(id);
        config.setName(displayName);


        if (genomeDescriptor.hasCytobands()) {
            InputStream cytobandStream = null;
            try {
                cytobandStream = genomeDescriptor.getCytoBandStream();
                BufferedReader reader = new BufferedReader(new InputStreamReader(cytobandStream));
                config.setCytobands(CytoBandFileParser.loadData(reader));

            } catch (IOException ex) {
                log.error("Error loading cytoband file", ex);
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
                config.setChromAliases(ChromAliasParser.loadChrAliases  (reader));
            }
        } catch (IOException e) {
            // We don't want to bomb if the alias load fails.  Just log it and proceed.
            log.error("Error loading chromosome alias table");
        } finally {
            closeSilently(aliasStream);
        }

        Genome   newGenome = new Genome(config);


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

}
