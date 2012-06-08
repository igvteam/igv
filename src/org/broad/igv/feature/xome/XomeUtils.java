package org.broad.igv.feature.xome;

import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.Exon;
import org.broad.igv.feature.FeatureUtils;
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
public class XomeUtils {


    static Map<String, List<Block>> blockCache = new HashMap();

    public static final Comparator<Block> GENOME_POSITION_COMPARATOR = new Comparator<Block>() {
        public int compare(Block o1, Block o2) {
            int genomeStart2 = o2.getGenomeStart();
            int genomeStart1 = o1.getGenomeStart();
            if (genomeStart2 >= genomeStart1 && o2.getGenomeEnd() <= o1.getGenomeEnd()) {
                return 0;
            } else {
                return genomeStart1 - genomeStart2;
            }
        }
    };

    public static final Comparator<Block> EXOME_POSITION_COMPARATOR = new Comparator<Block>() {
        public int compare(Block o1, Block o2) {
            int exomeStart2 = o2.getExomeStart();
            int exomeStart1 = o1.getExomeStart();
            if (exomeStart2 >= exomeStart1 && o2.getExomeEnd() <= o1.getExomeEnd()) {
                return 0;
            } else {
                return exomeStart1 - exomeStart2;
            }
        }
    };

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
        Block block = new Block(blockIndex, f.getStart(), exomeStart, f.getEnd() - f.getStart());

        while (iter.hasNext()) {
            f = iter.next();
            if (f.getStart() > block.getGenomeEnd()) {
                blocks.add(block);

                // Start next block
                int blockLength = block.getLength();
                exomeStart += blockLength;
                blockIndex++;
                int startingLength = f.getEnd() - f.getStart();
                block = new Block(blockIndex, f.getStart(), exomeStart, startingLength);
            } else {
                block.extend(f.getEnd());
            }

        }

        return blocks;

    }

    public static int genomeToExomePosition(List<Block> blocks, int genomePosition) {
        int idx = getIndexForGenomePosition(blocks, genomePosition);
        Block b = blocks.get(idx);

        if (genomePosition < b.getGenomeStart()) {
            return b.getExomeStart();
        } else {
            return genomePosition < b.getGenomeEnd() ?
                    b.getExomeStart() + (genomePosition - b.getGenomeStart()) :
                    b.getExomeEnd();
        }
    }


    public static int exomeToGenomePosition(List<Block> blocks, int exomePosition) {

        Block b = getBlockAtExomePosition(blocks, exomePosition);
        if (b != null) {
            return b.getGenomeStart() + (exomePosition - b.getExomeStart());
        } else {
            // ?? Should be impossible if "exomePosition" is in bounds
            b = blocks.get(blocks.size() - 1);
            return b.getGenomeEnd();
        }
    }

    public static Block getBlockAtGenomePosition(List<Block> blocks, int genomePosition) {
        Block key = new Block(-1, genomePosition, -1, 1);
        int r = Collections.binarySearch(blocks, key, GENOME_POSITION_COMPARATOR);
        if (r >= 0) {
            return blocks.get(r);
        } else {
            return null;
        }
    }


    public static Block getBlockAtExomePosition(List<Block> blocks, int exomePosition) {
        Block key = new Block(-1, -1, exomePosition, 1);
        int r = Collections.binarySearch(blocks, key, EXOME_POSITION_COMPARATOR);
        if (r >= 0) {
            return blocks.get(r);
        } else {
            return null;
        }
    }


    /**
     * Return the index to the last feature in the list with a start < the given position
     *
     * @param position
     * @param blocks
     * @return
     */
    public static int getIndexForGenomePosition(List<Block> blocks, double position) {


        int startIdx = 0;
        int endIdx = blocks.size() - 1;

        while (startIdx != endIdx) {
            int idx = (startIdx + endIdx) / 2;
            double distance = blocks.get(idx).getGenomeStart() - position;
            if (distance <= 0) {
                startIdx = idx;
            } else {
                endIdx = idx;
            }
            if (endIdx - startIdx < 10) {
                break;
            }
        }

        if (blocks.get(endIdx).getGenomeStart() >= position) {
            for (int idx = endIdx; idx >= 0; idx--) {
                if (blocks.get(idx).getGenomeStart() <= position) {
                    return idx;
                }
            }
            return 0;
        } else {
            for (int idx = endIdx + 1; idx < blocks.size(); idx++) {
                if (blocks.get(idx).getGenomeStart() >= position) {
                    return idx - 1;
                }

            }
            return blocks.size() - 1;
        }
    }


    public static int getIndexForExomePosition(List<Block> blocks, int exomeEnd) {
        return (getBlockAtExomePosition(blocks, exomeEnd).getIdx());
    }
}
