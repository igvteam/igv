package org.broad.igv.feature.xome;

import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.ui.panel.ReferenceFrame;

import java.util.List;

/**
 * @author Jim Robinson
 * @date 5/30/12
 */
public class ExomeReferenceFrame extends ReferenceFrame {


    // The screen pixel gap between blocks
    public static int blockGap = 2;

    List<Block> blocks;
    int firstBlockIdx;
    int endBlockIdx;

    int exomeOrigin;
    int genomeEnd;

    int exomeMaxPosition;

    public ExomeReferenceFrame(ReferenceFrame otherFrame, List<Block> blocks) {
        super(otherFrame);
        this.blocks = blocks;
        exomeMaxPosition = blocks.get(blocks.size() - 1).getExomeEnd();

        int idx = FeatureUtils.getIndexBefore(origin, blocks);
        if (blocks.get(idx).compareGenomePosition(origin) == 0) {
            firstBlockIdx = idx;
        } else {
            firstBlockIdx = (idx + 1) < blocks.size() ? (idx + 1) : idx;
        }

        findEnd();
    }


    @Override
    public void shiftOriginPixels(double delta) {

        if (exomeOrigin == 0 && delta < 0 || exomeOrigin >= exomeMaxPosition && delta > 0) return;

        double shiftBP = delta * getScale();
        exomeOrigin += shiftBP;
        if (exomeOrigin < 0) exomeOrigin = 0;
        if (exomeOrigin > exomeMaxPosition) exomeOrigin = exomeMaxPosition;

        // Find exome block that contains the new position.  We're assuming is very close to the current block.
        Block b = blocks.get(firstBlockIdx);
        int comp = b.compareExomePosition(exomeOrigin);
        if (comp > 0) {
            while (firstBlockIdx < blocks.size() - 1) {
                firstBlockIdx++;
                b = blocks.get(firstBlockIdx);
                if (b.compareExomePosition(exomeOrigin) <= 0) break;
            }
        } else if (comp < 0) {
            while (firstBlockIdx > 0) {
                firstBlockIdx--;
                b = blocks.get(firstBlockIdx);
                if (b.compareExomePosition(exomeOrigin) >= 0) break;
            }
        }

        // Find genomePosition
        double genomePosition = b.getGenomeStart() + (exomeOrigin - b.getExomeStart());

        super.setOrigin(genomePosition, true);
    }

    @Override
    public void snapToGrid() {
        super.setOrigin(Math.round(origin), true);
    }

    @Override
    public void setOrigin(double genomePosition, boolean repaint) {

        // Find the exomic block containing the genome position.  No assumption  made re close to current block.

        int idx = FeatureUtils.getIndexBefore(genomePosition, blocks);
        if (blocks.get(idx).compareGenomePosition(genomePosition) == 0) {
            firstBlockIdx = idx;
        } else {
            firstBlockIdx = (idx + 1) < blocks.size() ? (idx + 1) : idx;
        }
        findEnd();

        super.setOrigin(genomePosition, repaint);


    }

    private void findEnd() {
        double bpExtent = widthInPixels * getScale();

        int bp = 0;
        int idx = firstBlockIdx;
        Block firstBlock = blocks.get(firstBlockIdx);
        Block lastBlock = firstBlock;
        while (idx < blocks.size()) {
            lastBlock = blocks.get(idx);
            bp += lastBlock.getLength();
            if (bp > bpExtent) {
                break;
            }
            idx++;
        }


        exomeOrigin = firstBlock.getExomeStart() + (int) (origin - firstBlock.getGenomeStart());
        genomeEnd = lastBlock.getGenomeEnd();
    }

    /**
     * Return the chromosome (genomic) position corresponding to the screen pixel position.
     *
     * @param screenPosition
     * @return
     */
    @Override
    public double getChromosomePosition(int screenPosition) {

        int idx = firstBlockIdx;
        Block b;
        do {
            b = blocks.get(idx);
            int rightPixel = b.getRightPixel();
            if (rightPixel > screenPosition) {
                double delta = (screenPosition - b.getLeftPixel()) * getScale();
                return b.getGenomeStart() + delta;
            }
            idx++;
        } while (idx < blocks.size());
        return -1;

    }

    @Override
    public double getEnd() {
        return genomeEnd;
    }

    public List<Block> getBlocks() {
        return blocks;
    }

    public int getFirstBlockIdx() {
        return firstBlockIdx;
    }

    public boolean isExomeMode() {
        return true;
    }

    public int getExomeOrigin() {
        return exomeOrigin;
    }
}
