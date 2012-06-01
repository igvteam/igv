package org.broad.igv.feature.xome;

import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.Exon;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;

import java.io.IOException;
import java.util.*;

/**
 * @author Jim Robinson
 * @date 5/24/12
 */
public class XomeUtils {


    static Map<String, List<Block>> blockCache = new HashMap();

    public static synchronized List<Block> getBlocks(String chr) {

        List<Block> blocks = blockCache.get(chr);
        if (blocks == null) {
            Genome genome = GenomeManager.getInstance().getCurrentGenome();
            FeatureTrack geneTrack = IGV.getInstance().getGeneTrack();
            Chromosome chromosome = genome.getChromosome(chr);
            List<Feature> features = geneTrack.getFeatures(chr, 0, chromosome.getLength());
            blocks = collapseTranscripts(features);
            blockCache.put(chr, blocks);
        }
        return blocks;

    }


    public static List<Block> collapseTranscripts(List<Feature> features) {

        // Step 1,  extract exons
        List<Feature> exons = new ArrayList<Feature>(features.size() * 10);
        for (Feature f : features) {
            if (f instanceof BasicFeature) {
                List<Exon> tmp = ((BasicFeature) f).getExons();
                if (tmp != null) {
                    exons.addAll(tmp);
                }
            } else {
                // Dont know what to do with it
            }
        }

        FeatureUtils.sortFeatureList(exons);

        return collapseFeatures(exons);

    }

    /**
     * Collapse the list of features into a list of non-overlapping blocks.  It is assumed that the
     * feature list is sorted by start position.
     *
     * @param features
     * @return
     */
    public static List<Block> collapseFeatures(List<Feature> features) {

        List<Block> blocks = new ArrayList();
        if (features.isEmpty()) {
            return blocks;
        }

        Iterator<Feature> iter = features.iterator();
        Feature f = iter.next();
        int blockIndex = 0;
        int exomeStart = 0;
        Block block = new Block(blockIndex, f.getStart(), f.getEnd(), exomeStart);

        while (iter.hasNext()) {
            f = iter.next();
            if (f.getStart() > block.getGenomeEnd()) {
                blocks.add(block);

                // Start next block
                int blockLength = block.getLength();
                exomeStart += blockLength;
                blockIndex++;
                block = new Block(blockIndex, f.getStart(), f.getEnd(), exomeStart);
            } else {
                block.extend(f.getEnd());
            }

        }

        return blocks;

    }
}
