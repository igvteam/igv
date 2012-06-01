package org.broad.igv.feature.xome;

import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.ui.panel.ReferenceFrame;

import java.util.List;

/**
 * @author Jim Robinson
 * @date 5/30/12
 */
public class ExomeReferenceFrame extends ReferenceFrame {


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
        findEnd();
    }


    @Override
    public void shiftOriginPixels(double delta) {

        if (exomeOrigin == 0 || exomeOrigin >= exomeMaxPosition) return;

        double shiftBP = delta * getScale();
        exomeOrigin += shiftBP;
        if (exomeOrigin < 0) exomeOrigin = 0;
        if (exomeOrigin > exomeMaxPosition) exomeOrigin = exomeMaxPosition;

        // Find exome block that contains the new position
        Block b = blocks.get(firstBlockIdx);
        int comp = b.compareExomePosition(exomeOrigin);
        if (comp > 0) {
            while (firstBlockIdx < blocks.size() - 1) {
                firstBlockIdx++;
                b = blocks.get(firstBlockIdx);
                if (b.compareExomePosition(exomeOrigin) == 0) break;
            }
        } else if (comp < 0) {
            while (firstBlockIdx > 0) {
                firstBlockIdx--;
                b = blocks.get(firstBlockIdx);
                if (b.compareExomePosition(exomeOrigin) == 0) break;
            }
        }

        // Find genomePosition
        double genomePosition = b.getGenomeStart() + (exomeOrigin - b.getExomeStart());

        super.setOrigin(genomePosition);
    }

    @Override
    public void setOrigin(double genomePosition, boolean repaint) {

        // Find the exomic block containing the genome position
        int idx = FeatureUtils.getIndexBefore(genomePosition, blocks);
        if (blocks.get(idx).compareGenomePosition(genomePosition) == 0) {
            firstBlockIdx = idx;
        } else {
            firstBlockIdx = (idx + 1) < blocks.size() ? (idx + 1) : idx;
        }

        // double adjustedOrigin =

        super.setOrigin(genomePosition, repaint);

        findEnd();
    }

    private void findEnd() {
        double bpExtent = widthInPixels * getScale();
        Block firstBlock = (Block) FeatureUtils.getFeatureAfter(origin, blocks);
        Block lastBlock = firstBlock;

        int bp = 0;
        int idx = firstBlock.getIdx();
        while (idx < blocks.size()) {
            lastBlock = blocks.get(idx);
            bp += lastBlock.getLength();
            if (bp > bpExtent) {
                break;
            }
            idx++;
        }
        exomeOrigin = firstBlock.getExomeStart() +  (int) (origin - firstBlock.getGenomeStart());
        genomeEnd = lastBlock.getGenomeEnd();
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
}
