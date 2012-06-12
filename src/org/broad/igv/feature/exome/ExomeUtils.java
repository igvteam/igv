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


    static Map<String, ExomeData> blockCache = new HashMap();


    private static final Comparator<ExomeBlock> GENOME_POSITION_COMPARATOR = new Comparator<ExomeBlock>() {
        public int compare(ExomeBlock o1, ExomeBlock o2) {
            int genomeStart2 = o2.getGenomeStart();
            int genomeStart1 = o1.getGenomeStart();
            if (genomeStart2 >= genomeStart1 && o2.getGenomeEnd() <= o1.getGenomeEnd()) {
                return 0;
            } else {
                return genomeStart1 - genomeStart2;
            }
        }
    };

    private static final Comparator<ExomeBlock> EXOME_POSITION_COMPARATOR = new Comparator<ExomeBlock>() {
        public int compare(ExomeBlock o1, ExomeBlock o2) {
            int exomeStart2 = o2.getExomeStart();
            int exomeStart1 = o1.getExomeStart();
            if (exomeStart2 >= exomeStart1 && o2.getExomeEnd() <= o1.getExomeEnd()) {
                return 0;
            } else {
                return exomeStart1 - exomeStart2;
            }
        }
    };

    public static List<ExomeBlock> getBlocks(String chr) {
        ExomeData exomeData = getExomeData(chr);
        return exomeData.block;
    }

    public static List<Gene> getGenes(String chr) {
        ExomeData exomeData = getExomeData(chr);
        return exomeData.genes;

    }

    private static synchronized ExomeData getExomeData(String chr) {
        ExomeData exomeData = blockCache.get(chr);
        if (exomeData == null) {
            Genome genome = GenomeManager.getInstance().getCurrentGenome();
            FeatureTrack geneTrack = IGV.getInstance().getGeneTrack();
            Chromosome chromosome = genome.getChromosome(chr);

            List<Feature> features = geneTrack.getFeatures(chr, 0, chromosome.getLength());
            List<ExomeBlock> blocks = collapseTranscripts(features);
            List<Gene> genes = collapseToGenes(features);
            exomeData = new ExomeData(blocks, genes);
            blockCache.put(chr, exomeData);

        }
        return exomeData;
    }


    private static List<Gene> collapseToGenes(List<Feature> features) {

        Map<String, Gene> genes = new HashMap<String, Gene>(10000);

        for (Feature f : features) {
            if (f instanceof BasicFeature) {
                String geneName = ((BasicFeature) f).getName();
                Gene gene = genes.get(geneName);
                if (gene == null) {
                    gene = new Gene((BasicFeature) f);
                    genes.put(geneName, gene);
                } else {
                    gene.expand(f);
                }
            } else {

            }
        }

        List<Gene> geneList = new ArrayList<Gene>(genes.values());
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

                // Start next block
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

    public static int genomeToExomePosition(List<ExomeBlock> blocks, int genomePosition) {
        int idx = getIndexForGenomePosition(blocks, genomePosition);
        ExomeBlock b = blocks.get(idx);

        if (genomePosition < b.getGenomeStart()) {
            return b.getExomeStart();
        } else {
            return genomePosition < b.getGenomeEnd() ?
                    b.getExomeStart() + (genomePosition - b.getGenomeStart()) :
                    b.getExomeEnd();
        }
    }


    public static int exomeToGenomePosition(List<ExomeBlock> blocks, int exomePosition) {

        ExomeBlock b = getBlockAtExomePosition(blocks, exomePosition);
        if (b != null) {
            return b.getGenomeStart() + (exomePosition - b.getExomeStart());
        } else {
            // ?? Should be impossible if "exomePosition" is in bounds
            b = blocks.get(blocks.size() - 1);
            return b.getGenomeEnd();
        }
    }

    private static ExomeBlock getBlockAtGenomePosition(List<ExomeBlock> blocks, int genomePosition) {
        ExomeBlock key = new ExomeBlock(-1, genomePosition, -1, 1);
        int r = Collections.binarySearch(blocks, key, GENOME_POSITION_COMPARATOR);
        if (r >= 0) {
            return blocks.get(r);
        } else {
            return null;
        }
    }


    private static ExomeBlock getBlockAtExomePosition(List<ExomeBlock> blocks, int exomePosition) {
        ExomeBlock key = new ExomeBlock(-1, -1, exomePosition, 1);
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
    public static int getIndexForGenomePosition(List<ExomeBlock> blocks, double position) {


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

    public static class Gene implements NamedFeature {

        String name;
        String chr;
        int start;
        int end;

        Gene(NamedFeature f) {
            name = f.getName();
            chr = f.getChr();
            start = f.getStart();
            end = f.getEnd();
        }

        @Override
        public String getChr() {
            return chr;  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        public int getStart() {
            return start;  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        public int getEnd() {
            return end;  //To change body of implemented methods use File | Settings | File Templates.
        }


        @Override
        public String getName() {
            return name;
        }

        public void expand(Feature f) {
            //if(!chr.equals(f.getChr())) throw error
            start = Math.min(start, f.getStart());
            end = Math.max(end, f.getEnd());
        }

    }

    static class ExomeData {
        List<ExomeBlock> block;
        List<Gene> genes;

        ExomeData(List<ExomeBlock> block, List<Gene> genes) {
            this.block = block;
            this.genes = genes;
        }

        public List<ExomeBlock> getBlock() {
            return block;
        }

        public List<Gene> getGenes() {
            return genes;
        }
    }
}
