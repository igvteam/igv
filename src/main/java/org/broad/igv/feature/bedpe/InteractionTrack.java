package org.broad.igv.feature.bedpe;

import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.sprite.ClusterTrack;
import org.broad.igv.track.AbstractTrack;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.TrackClickEvent;
import org.broad.igv.track.TrackMenuUtils;
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
 * Created by jrobinso on 6/29/18.
 */
public class InteractionTrack extends AbstractTrack {

    private Map<String, List<BedPEFeature>> featureMap;

    private PEArcRenderer renderer;

    public InteractionTrack(ResourceLocator locator, List<BedPEFeature> featureList, Genome genome) {
        super(locator);
        init(featureList, genome);
        renderer = new PEArcRenderer();
        setHeight(250, true);
        setColor(new Color(180, 25, 137));
    }

    private void init(List<BedPEFeature> featureList, Genome genome) {

        this.featureMap = new HashMap<>();

        // Sort feature lists by "start" (minimum of start1, start2)
        Collections.sort(featureList, (o1, o2) -> o1.getStart() - o2.getStart());

        for (BedPEFeature f : featureList) {

            String key;
            if (f.chr1.equals(f.chr2)) {
                key = genome == null ? f.chr1 : genome.getCanonicalChrName(f.chr1);
            } else {
                key = "OTHER";
            }

            List<BedPEFeature> features = featureMap.get(key);
            if (features == null) {
                features = new ArrayList<>();
                featureMap.put(key, features);
            }

            features.add(f);
        }

        if(featureMap.containsKey("OTHER")) {
            featureMap.put(Globals.CHR_ALL, createWGFeatures(featureMap.get("OTHER"), genome));
        }
    }

    private List<BedPEFeature> createWGFeatures(List<BedPEFeature> features, Genome genome) {

        List<BedPEFeature> wgFeatures = new ArrayList<>(features.size());

        for (BedPEFeature f : features) {

            BedPEFeature wgFeature = new BedPEFeature();
            wgFeature.chr1 = Globals.CHR_ALL;
            wgFeature.chr2 = Globals.CHR_ALL;
            wgFeature.name = f.name;
            wgFeature.score = f.score;
            wgFeature.start1 = genome.getGenomeCoordinate(f.chr1, f.start1);
            wgFeature.end1 = genome.getGenomeCoordinate(f.chr1, f.end1);
            wgFeature.start2 = genome.getGenomeCoordinate(f.chr1, f.start2);
            wgFeature.end2 = genome.getGenomeCoordinate(f.chr1, f.end2);
            wgFeatures.add(wgFeature);

        }
        return wgFeatures;
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

        Graphics2D g2d = context.getGraphics();
        Rectangle clip = new Rectangle(g2d.getClip().getBounds());
        g2d.setClip(rect.intersection(clip.getBounds()));
        context.clearGraphicsCache();

        try {
            String chr = context.getReferenceFrame().getChrName();
            List<BedPEFeature> features = featureMap.get(chr);

            if (features != null) {
                renderer.render(features, context, rect, this);
            }
            context.clearGraphicsCache();
        } finally {
            g2d.setClip(clip);
        }
    }


    @Override
    public IGVPopupMenu getPopupMenu(TrackClickEvent te) {

        IGVPopupMenu menu = new IGVPopupMenu();

        menu.add(TrackMenuUtils.getTrackRenameItem(Collections.singleton(InteractionTrack.this)));

        JMenuItem item = new JMenuItem("Set Track Height...");
        item.addActionListener(evt -> TrackMenuUtils.changeTrackHeight(Collections.singleton(InteractionTrack.this)));
        menu.add(item);

        item = new JMenuItem("Set Track Color...");
        item.addActionListener(evt -> TrackMenuUtils.changeTrackColor(Collections.singleton(InteractionTrack.this)));
        menu.add(item);

        item = new JMenuItem("Toggle Arc Direction");
        item.addActionListener(evt -> {
            if (renderer.direction == PEArcRenderer.Direction.UP) {
                renderer.direction = PEArcRenderer.Direction.DOWN;
            } else {
                renderer.direction = PEArcRenderer.Direction.UP;
            }
            IGV.getInstance().repaint();
        });
        menu.add(item);

        item = new JMenuItem("Set Line Thickness...");
        item.addActionListener(e -> {
            String t = MessageUtils.showInputDialog("Line thickness", String.valueOf(renderer.lineThickness));
            if (t != null) {
                try {
                    renderer.lineThickness = Integer.parseInt(t);
                    IGV.getInstance().repaint();
                } catch (NumberFormatException e1) {
                    MessageUtils.showErrorMessage("Line thickness must be an integer", e1);
                }
            }
        });
        menu.add(item);

        return menu;
    }
}
