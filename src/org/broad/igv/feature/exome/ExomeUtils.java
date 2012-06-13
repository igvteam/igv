package org.broad.igv.feature.exome;

import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.ui.IGV;
import org.broad.tribble.Feature;

import java.util.*;

/**
 * @author Jim Robinson
 * @date 5/24/12
 */
public class ExomeUtils {


    public static List<ExomeReferenceFrame.Gene> collapseToGenes(List<Feature> features) {

        Map<String, ExomeReferenceFrame.Gene> genes = new HashMap<String, ExomeReferenceFrame.Gene>(10000);

        for (Feature f : features) {
            if (f instanceof BasicFeature) {
                String geneName = ((BasicFeature) f).getName();
                ExomeReferenceFrame.Gene gene = genes.get(geneName);
                if (gene == null) {
                    gene = new ExomeReferenceFrame.Gene((BasicFeature) f);
                    genes.put(geneName, gene);
                } else {
                    gene.expand(f);
                }
            } else {

            }
        }

        List<ExomeReferenceFrame.Gene> geneList = new ArrayList<ExomeReferenceFrame.Gene>(genes.values());
        FeatureUtils.sortFeatureList(geneList);
        return geneList;
    }

    static List<ExomeBlock> collapseTranscripts(List<Feature> features) {

        // Step 1,  extract exons
        List<Feature> exons = new ArrayList<Feature>(features.size() * 10);
        for (Feature f : features) {
            if (f instanceof BasicFeature) {
                List<Exon> tmp = ((BasicFeature) f).getExons();
                if (tmp != null) {
                    exons.addAll(tmp);
                }
            } else {

            }
        }

        FeatureUtils.sortFeatureList(exons);
        List<ExomeBlock> blocks = collapseFeatures(exons);
        return blocks;

    }

    /**
     * Collapse the list of features into a list of non-overlapping blocks.  It is assumed that the
     * feature list is sorted by start position.
     *
     * @param features
     * @return
     */
    private static List<ExomeBlock> collapseFeatures(List<Feature> features) {

        List<ExomeBlock> blocks = new ArrayList();
        if (features.isEmpty()) {
            return blocks;
        }

        Iterator<Feature> iter = features.iterator();
        Feature f = iter.next();
        int blockIndex = 0;
        int exomeStart = 0;
        ExomeBlock block = new ExomeBlock(blockIndex, f.getStart(), exomeStart, f.getEnd() - f.getStart());

        while (iter.hasNext()) {
            f = iter.next();
            if (f.getStart() > block.getGenomeEnd()) {
                blocks.add(block);

                // Start next blocks
                int blockLength = block.getLength();
                exomeStart += blockLength;
                blockIndex++;
                int startingLength = f.getEnd() - f.getStart();
                block = new ExomeBlock(blockIndex, f.getStart(), exomeStart, startingLength);
            } else {
                block.extend(f.getEnd());
            }

        }

        return blocks;

    }
}
