package org.broad.igv.feature.xome;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.feature.FeatureUtils;
import org.broad.igv.feature.Locus;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.DragEventManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.MessageUtils;

import java.util.List;

/**
 * @author Jim Robinson
 * @date 5/30/12
 */
public class ExomeReferenceFrame extends ReferenceFrame {

    private static Logger log = Logger.getLogger(ExomeReferenceFrame.class);

    // The screen pixel gap between blocks
    // public static int blockGap = 2;

    //List<Block> blocks;
    int firstBlockIdx;
    int endBlockIdx;

    int exomeOrigin;
    int exomeEnd;

    int genomeEnd;
    private int blockGap;

    public ExomeReferenceFrame(ReferenceFrame otherFrame) {

        super(otherFrame);

    }

    private void findFirstBlock(List<Block> blocks) {
        int idx = FeatureUtils.getIndexBefore(origin, blocks);
        final int comp = blocks.get(idx).compareGenomePosition(origin);
        if (comp == 0) {
            firstBlockIdx = idx;
        } else {
            firstBlockIdx = (idx + 1) < blocks.size() ? (idx + 1) : idx;
        }
    }


    @Override
    public void shiftOriginPixels(double delta) {

        if (exomeOrigin == 0 && delta < 0) return;

        double shiftBP = delta * getScale();
        exomeOrigin += shiftBP;
        if (exomeOrigin < 0) exomeOrigin = 0;

        // Find exome block that contains the new position.  We're assuming is very close to the current block.
        List<Block> blocks = XomeUtils.getBlocks(getChrName());
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
        List<Block> blocks = XomeUtils.getBlocks(getChrName());
        int idx = FeatureUtils.getIndexBefore(genomePosition, blocks);
        if (blocks.get(idx).compareGenomePosition(genomePosition) == 0) {
            firstBlockIdx = idx;
        } else {
            firstBlockIdx = (idx + 1) < blocks.size() ? (idx + 1) : idx;
        }
        findEnd(blocks);

        super.setOrigin(genomePosition, repaint);


    }

    @Override
    public void jumpTo(String chr, int start, int end) {

        setInterval(new Locus(chr, start, end));
        IGV.getInstance().repaintDataAndHeaderPanels();
        IGV.getInstance().repaintStatusAndZoomSlider();

    }

    /**
     * Jump to a specific locus (in genome coordinates).
     *
     * @param locus
     */
    @Override
    public void setInterval(Locus locus) {

        this.initialLocus = locus;
        this.chrName = locus.getChr();
        this.origin = locus.getStart();


        List<Block> blocks = XomeUtils.getBlocks(chrName);
        findFirstBlock(blocks);
        Block firstBlock = blocks.get(firstBlockIdx);

        origin = firstBlock.getGenomeStart();
        exomeOrigin = firstBlock.getExomeStart();

        genomeEnd = locus.getEnd();
        // Loop through blocks until we find the genome end

        Block endBlock;
        int idx = firstBlockIdx;
        while (idx < blocks.size()) {
            endBlock = blocks.get(idx);
            if (endBlock.getGenomeStart() > genomeEnd) {
                endBlockIdx = idx - 1;
                break;
            } else if (endBlock.getGenomeEnd() >= genomeEnd) {

                endBlockIdx = idx;
                break;
            }
            idx++;
        }
        genomeEnd = blocks.get(endBlockIdx).getGenomeEnd();

        // Get total base pairs traversed, in exome coordinates, and # of pixels available
        int nBlocks = (endBlockIdx - firstBlockIdx) + 1;
        blockGap = nBlocks < 80 ? 2 : (nBlocks < 160 ? 1 : 0);

        int pWidth = widthInPixels > 0 ? widthInPixels : 1000; // <= if not known guess
        pWidth -= 30;
        int bp = 0;
        for (idx = firstBlockIdx + 1; idx <= endBlockIdx; idx++) {
            bp += blocks.get(idx).getLength();
            pWidth -= blockGap;
        }

        if (pWidth < 10) pWidth = 10;

        locationScale = ((double) bp) / pWidth;
        locationScaleValid = true;

        exomeEnd = exomeOrigin + bp;

        imputeZoom(exomeOrigin, exomeEnd);
    }

    @Override
    public void zoomTo(int newZoom, double newCenter) {

        // Zoom in exome coordinates.  This is approximate, we aren't accounting for gaps
        int exomeLength = exomeEnd - exomeOrigin;
        exomeOrigin += exomeLength / 4;
        exomeEnd -= exomeLength / 4;
        locationScale /= 2;
        locationScaleValid = true;

        List<Block> blocks = XomeUtils.getBlocks(chrName);
        for (int idx = firstBlockIdx; idx < blocks.size(); idx++) {
            if (blocks.get(idx).getExomeEnd() > exomeOrigin) {
                firstBlockIdx = idx;
                break;
            }
        }
        Block b = blocks.get(firstBlockIdx);
        origin = b.getGenomeStart() + (exomeOrigin - b.getExomeStart());

        for (int idx = endBlockIdx; idx > firstBlockIdx; idx--) {
            if (blocks.get(idx).getExomeStart() < exomeEnd) {
                endBlockIdx = idx;
                break;
            }
        }
        b = blocks.get(endBlockIdx);
        genomeEnd = b.getGenomeEnd() + (exomeEnd - b.getExomeStart());

        int nBlocks = (endBlockIdx - firstBlockIdx) + 1;
        blockGap = nBlocks < 80 ? 2 : (nBlocks < 160 ? 1 : 0);

        imputeZoom(exomeOrigin, exomeEnd);

        recordHistory();

        IGV.getInstance().repaintDataAndHeaderPanels();
        IGV.getInstance().repaintStatusAndZoomSlider();

        // This is a hack,  this is not a drag event but is a "jump"
        DragEventManager.getInstance().dragStopped();
    }

    private void findEnd(List<Block> blocks) {
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

        exomeEnd = lastBlock.getExomeEnd();

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

        List<Block> blocks = XomeUtils.getBlocks(getChrName());
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
        return XomeUtils.getBlocks(getChrName());
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

    public int getBlockGap() {
        return blockGap;
    }
}
