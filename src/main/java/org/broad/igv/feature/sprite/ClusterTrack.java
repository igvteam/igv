package org.broad.igv.feature.sprite;

import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.*;
import java.util.List;

/**
 * Created by jrobinso on 6/30/18.
 */
public class ClusterTrack extends AbstractTrack {


    List<Cluster> clusters;
    int rowHeight = 2;


    public ClusterTrack(ResourceLocator locator, List<Cluster> clusters, Genome genome) {
        super(locator);
        this.clusters = clusters;
        computeWGScores(clusters, genome);
    }

    @Override
    public int getHeight() {
        return clusters.size() * rowHeight;
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

        int y = 0;

        for (Cluster c : clusters) {

            List<Integer> loci = c.posMap.get(chr);

            if (loci != null) {
                for (Integer position : loci) {

                    double pixelStart = ((position - origin) / locScale);
                    double pixelEnd = ((position - origin) / locScale);

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

                    y += rowHeight;
                }
            }
        }
    }

    private void computeWGScores(List<Cluster> clusters, Genome genome) {

        // Bin whole genome
        int nBins = 1000;
        double binSize = (genome.getTotalLength() / 1000.0) / nBins;

        for (Cluster cluster : clusters) {

            int[] bins = new int[nBins];
            List<Integer> occupiedBins = new ArrayList<>();

            for (Map.Entry<String, List<Integer>> entry : cluster.posMap.entrySet()) {

                String chr = entry.getKey();
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

        final JMenuItem rowHeightItem = new JMenuItem("Set Row Height...");
        rowHeightItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                String t = MessageUtils.showInputDialog("Row height", String.valueOf(rowHeight));
                if(t != null) {
                    try {
                        int h = Integer.parseInt(t);
                        rowHeight = h;
                        IGV.getInstance().repaint();
                    } catch (NumberFormatException e1) {
                        MessageUtils.showErrorMessage("Row height must be a number", e1);
                    }
                }
            }
        });
        menu.add(rowHeightItem);

        JMenuItem item = new JMenuItem("Set Track Color...");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                TrackMenuUtils.changeTrackColor(Collections.singleton(ClusterTrack.this));
            }
        });
        menu.add(item);

        return menu;
    }

    @Override
    public String getValueStringAt(String chr, double position, int mouseX, int mouseY, ReferenceFrame frame) {

        int row = mouseY / rowHeight;

        if(row < clusters.size()) {
            return clusters.get(row).name;
        }
        else {
            return "";
        }

    }
}