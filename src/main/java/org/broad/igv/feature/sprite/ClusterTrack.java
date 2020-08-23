package org.broad.igv.feature.sprite;

import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.AbstractTrack;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.TrackClickEvent;
import org.broad.igv.track.TrackMenuUtils;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ResourceLocator;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import javax.swing.*;
import java.awt.*;
import java.util.*;
import java.util.List;

/**
 * Created by jrobinso on 6/30/18.
 */
public class ClusterTrack extends AbstractTrack {


    int binSize;
    ClusterParser.ClusterSet clusterSet;
    List<Cluster> binnedClusters;
    List<Cluster> wgClusters;
    int rowHeight = 2;
    Genome genome;

    public ClusterTrack() {
    }

    public ClusterTrack(ResourceLocator locator, ClusterParser.ClusterSet clusterSet, Genome genome) {
        super(locator);
        this.clusterSet = clusterSet;
        this.binSize = clusterSet.binSize;
        this.genome = genome;
        this.binnedClusters = computeBinnedClusters(clusterSet.binSize, clusterSet.clusters, genome);
        computeWGScores(this.binnedClusters, genome);
    }

    private List<Cluster> computeBinnedClusters(int binSize, List<Cluster> clusters, Genome genome) {

        List<Cluster> binnedClusters = new ArrayList<>();

        for (Cluster c : clusters) {

            String name = c.name;
            Map<String, List<Integer>> posMap = new HashMap<>();
            for (Map.Entry<String, List<Integer>> entry : c.posMap.entrySet()) {

                String chr = genome.getCanonicalChrName(entry.getKey());
                Set<Integer> bins = new HashSet<>();
                for (Integer pos : entry.getValue()) {
                    bins.add((pos / binSize) * binSize);
                }
                List<Integer> binList = new ArrayList<>(bins);
                Collections.sort(binList);
                posMap.put(chr, binList);
            }
            binnedClusters.add(new Cluster(name, posMap));
        }

        return binnedClusters;
    }


    public void setBinSize(int binSize) {
        this.binSize = binSize;
        this.binnedClusters = computeBinnedClusters(this.binSize, clusterSet.clusters, this.genome);
        computeWGScores(this.binnedClusters, genome);
    }


    @Override
    public int getHeight() {
        return binnedClusters.size() * rowHeight;
    }

    @Override
    public boolean isReadyToPaint(ReferenceFrame frame) {
        return true;
    }

    @Override
    public void load(ReferenceFrame frame) {
        // Nothing to do, this track is pre-loaded
    }

    @Override
    public void render(RenderContext context, Rectangle rect) {

        String chr = context.getReferenceFrame().getChrName();
        double origin = context.getOrigin();
        double locScale = context.getScale();

        double binSize = chr.equals(Globals.CHR_ALL) ? this.binSize / 1000 : this.binSize;

        int y = 0;

        for (Cluster c : binnedClusters) {

            List<Integer> loci = c.posMap.get(chr);

            if (loci != null) {
                for (Integer position : loci) {

                    double pixelStart = ((position - origin) / locScale);
                    double pixelEnd = ((position + binSize - origin) / locScale);

                    // If the any part of the feature fits in the Track rectangle draw it
                    if (pixelEnd >= rect.getX() && pixelStart <= rect.getMaxX()) {

                        Color color = this.getColor();
                        Graphics2D g = context.getGraphic2DForColor(color);

                        int w = (int) (pixelEnd - pixelStart);
                        if (w < 3) {
                            w = 3;
                            pixelStart--;
                        }

                        g.fillRect((int) pixelStart, y, w, rowHeight);

                    }
                }
            }
            y += rowHeight;
        }
    }


    private void computeWGScores(List<Cluster> clusters, Genome genome) {

        // Bin whole genome
        int nBins = 1000;
        double binSize = (genome.getTotalLength() / 1000.0) / nBins;
        Set<String> wgChrNames = new HashSet<>(genome.getLongChromosomeNames());

        for (Cluster cluster : clusters) {

            int[] bins = new int[nBins];
            List<Integer> occupiedBins = new ArrayList<>();

            for (Map.Entry<String, List<Integer>> entry : cluster.posMap.entrySet()) {

                String chr = entry.getKey();

                if (!wgChrNames.contains(chr)) continue;

                List<Integer> posList = entry.getValue();

                for (Integer pos : posList) {

                    int genomeCoordinate = genome.getGenomeCoordinate(chr, pos);
                    int b = (int) (genomeCoordinate / binSize);
                    bins[b]++;
                }

            }

            for (int i = 0; i < bins.length; i++) {
                if (bins[i] > 0) {
                    occupiedBins.add((int) (i * binSize));
                }
            }
            cluster.posMap.put(Globals.CHR_ALL, occupiedBins);

        }

    }


    @Override
    public IGVPopupMenu getPopupMenu(TrackClickEvent te) {

        IGVPopupMenu menu = new IGVPopupMenu();

        menu.add(TrackMenuUtils.getTrackRenameItem(Collections.singleton(ClusterTrack.this)));

        final JMenuItem binSizeItem = new JMenuItem("Set Bin Size...");
        binSizeItem.addActionListener(e -> {
            String t = MessageUtils.showInputDialog("Bin Size", String.valueOf(ClusterTrack.this.binSize));
            if (t != null) {
                try {
                    int bs = Integer.parseInt(t);
                    ClusterTrack.this.setBinSize(bs);
                    IGV.getInstance().repaint();
                } catch (NumberFormatException e1) {
                    MessageUtils.showErrorMessage("Bin size must be an integer", e1);
                }
            }
        });
        menu.add(binSizeItem);

        final JMenuItem rowHeightItem = new JMenuItem("Set Row Height...");
        rowHeightItem.addActionListener(e -> {
            String t = MessageUtils.showInputDialog("Row height", String.valueOf(rowHeight));
            if (t != null) {
                try {
                    int h = Integer.parseInt(t);
                    rowHeight = h;
                    IGV.getInstance().repaint();
                } catch (NumberFormatException e1) {
                    MessageUtils.showErrorMessage("Row height must be a number", e1);
                }
            }
        });
        menu.add(rowHeightItem);

        JMenuItem item = new JMenuItem("Set Track Color...");
        item.addActionListener(evt -> TrackMenuUtils.changeTrackColor(Collections.singleton(ClusterTrack.this)));
        menu.add(item);

        return menu;
    }

    @Override
    public String getValueStringAt(String chr, double position, int mouseX, int mouseY, ReferenceFrame frame) {

        int row = mouseY / rowHeight;

        if (row < binnedClusters.size()) {
            return binnedClusters.get(row).name;
        } else {
            return "";
        }

    }


    @Override
    public void marshalXML(Document document, Element element) {

        super.marshalXML(document, element);

        element.setAttribute("binSize", String.valueOf(binSize));
        element.setAttribute("sequence", String.valueOf(rowHeight));

    }

    @Override
    public void unmarshalXML(Element element, Integer version) {

        super.unmarshalXML(element, version);

        this.binSize = Integer.parseInt(element.getAttribute("binSize"));
        this.rowHeight = Integer.parseInt(element.getAttribute("rowHeight"));

    }
}