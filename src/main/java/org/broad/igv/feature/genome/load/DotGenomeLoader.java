package org.broad.igv.feature.genome.load;

import htsjdk.tribble.Feature;
import org.apache.log4j.Logger;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.Sequence;
import org.broad.igv.feature.genome.SequenceWrapper;
import org.broad.igv.feature.genome.fasta.FastaBlockCompressedSequence;
import org.broad.igv.feature.genome.fasta.FastaIndexedSequence;
import org.broad.igv.feature.gff.GFFFeatureSource;
import org.broad.igv.track.FeatureCollectionSource;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.TrackType;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.io.*;
import java.util.Collection;
import java.util.List;

public class DotGenomeLoader extends GenomeLoader {

    private static Logger log = Logger.getLogger(DotGenomeLoader.class);

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
                geneFeatureTrack.setUrl(annotationURL);
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

        Genome newGenome;

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

}
