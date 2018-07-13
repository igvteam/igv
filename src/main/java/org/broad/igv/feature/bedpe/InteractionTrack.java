package org.broad.igv.feature.bedpe;

import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import javax.xml.bind.annotation.XmlAttribute;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.Arc2D;
import java.util.*;
import java.util.List;

import static org.broad.igv.feature.bedpe.InteractionTrack.Direction.DOWN;
import static org.broad.igv.feature.bedpe.InteractionTrack.Direction.UP;


/**
 * Created by jrobinso on 6/29/18.
 */
public class InteractionTrack extends AbstractTrack {

    enum Direction {UP, DOWN}

    @XmlAttribute
    InteractionTrack.Direction direction = DOWN;

    @XmlAttribute
    int thickness = 1;

    @XmlAttribute
    boolean  hideLargeFeatures = false;

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

        if (featureMap.containsKey("OTHER")) {
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
            if (direction == UP) {
                direction = Direction.DOWN;
            } else {
                direction = UP;
            }
            IGV.getInstance().repaint();
        });
        menu.add(item);

        item = new JMenuItem("Set Line Thickness...");
        item.addActionListener(e -> {
            String t = MessageUtils.showInputDialog("Line thickness", String.valueOf(thickness));
            if (t != null) {
                try {
                    thickness = Integer.parseInt(t);
                    IGV.getInstance().repaint();
                } catch (NumberFormatException e1) {
                    MessageUtils.showErrorMessage("Line thickness must be an integer", e1);
                }
            }
        });
        menu.add(item);

//        final JCheckBoxMenuItem cbItem = new JCheckBoxMenuItem("Hide Large Features");
//        cbItem.setSelected(hideLargeFeatures);
//        cbItem.addActionListener(e -> {
//            InteractionTrack.this.hideLargeFeatures = cbItem.isSelected();
//            IGV.getInstance().repaint();
//        });
//

       return menu;
    }


    /**
     * Created by jrobinso on 6/29/18.
     */
    public class PEArcRenderer {


        private Map<Color, Color> alphaColors = new HashMap<>();

        double theta = Math.toRadians(45);
        double sinTheta = Math.sin(theta);
        double cosTheta = Math.cos(theta);


        public void render(List<BedPEFeature> featureList, RenderContext context, Rectangle trackRectangle, Track track) {

            double origin = context.getOrigin();
            double locScale = context.getScale();

            // Autoscale theta
            double max = 0;
            for (BedPEFeature feature : featureList) {

                // Note -- don't cast these to an int until the range is checked.
                // could get an overflow.
                double pixelStart = ((feature.getStart() - origin) / locScale);
                double pixelEnd = ((feature.getEnd() - origin) / locScale);
                if (pixelEnd >= trackRectangle.getX() && pixelStart <= trackRectangle.getMaxX()) {
                    max = Math.max(max, pixelEnd - pixelStart);
                }
            }
            double a = Math.min(trackRectangle.width, max) / 2;
            if (max > 0) {
                double coa = trackRectangle.height / a;
                theta = SagittusEstimate.estimateTheta(coa);
                sinTheta = Math.sin(theta);
                cosTheta = Math.cos(theta);
            }


            Graphics2D g = (Graphics2D) context.getGraphics().create();
            Color trackColor = track.getColor();

            try {
                for (BedPEFeature feature : featureList) {

                    // Note -- don't cast these to an int until the range is checked.
                    // could get an overflow.
                    double pixelStart = ((feature.getStart() - origin) / locScale);
                    double pixelEnd = ((feature.getEnd() - origin) / locScale);
                    double width = pixelEnd - pixelStart;


                    // If the any part of the feature fits in the Track rectangle draw it
                    if (pixelEnd >= trackRectangle.getX() && pixelStart <= trackRectangle.getMaxX()) {

                        int s = (feature.start1 + feature.end1) / 2;
                        int e = (feature.start2 + feature.end2) / 2;

                        double ps = (s - origin) / locScale;
                        double pe = (e - origin) / locScale;

                        Color fcolor = feature.color == null ? trackColor : feature.color;
                        if(fcolor != null && width > trackRectangle.width) {
                            fcolor = getAlphaColor(fcolor);
                        }
                        if (fcolor != null) {
                            g.setColor(fcolor);
                        }
                        if (feature.thickness > 1) {
                            g.setStroke(new BasicStroke(feature.thickness));
                        }

                        drawArc(g, trackRectangle, track, ps, pe);

                    }
                }
            } finally {
                g.dispose();
            }
        }

        private Color getAlphaColor(Color fcolor) {
            Color ac = alphaColors.get(fcolor);
            if(ac == null) {
                float [] rgb = new float[3];
                rgb = fcolor.getRGBColorComponents(rgb);
                ac = new Color(rgb[0], rgb[1], rgb[2], 0.05f);
                alphaColors.put(fcolor, ac);
            }
            return ac;
        }

        private void drawArc(Graphics2D g, Rectangle trackRectangle, Track track, double x1, double x2) {

            double pixelStart = Math.min(x1, x2);
            double pixelEnd = Math.max(x1, x2);

            if (thickness > 1) g.setStroke(new BasicStroke(thickness));

            int w = (int) (pixelEnd - pixelStart);
            if (w < 3) {
                w = 3;
                pixelStart--;
            }

            double a = w / 2;
            double r = a / sinTheta;
            double b = cosTheta * r;
            //double c = r - b;       // sagittus
            //double d = r - a;       // x delta to corner of bounding box
            double x = pixelStart + a;
            final int trackBaseLine = trackRectangle.y + trackRectangle.height;
            double y = direction == UP ? trackBaseLine + b : trackRectangle.y - b;
            double angleSt = direction == UP ? 90 - Math.toDegrees(theta) : 270 - Math.toDegrees(theta);
            double ext = Math.toDegrees(2 * theta);

            Arc2D.Double arcPath = new Arc2D.Double();
            arcPath.setArcByCenter(x, y, r, angleSt, ext, Arc2D.OPEN);

            g.draw(arcPath);
        }


    }


    /**
     * Estimate theta given the ratio of track height to 1/2 the feature width (coa).  This relationship is approximately linear.
     */
    static class SagittusEstimate {

        static final double[] coa = {0.01570925532366355, 0.15838444032453644, 0.3249196962329063, 0.5095254494944288, 0.7265425280053609, 0.9999999999999999};
        static final double[] theta = {0.031415926535897934, 0.3141592653589793, 0.6283185307179586, 0.9424777960769379, 1.2566370614359172, 1.5707963267948966};

        static double estimateTheta(double x) {

            int idx;
            for (idx = 0; idx < coa.length; idx++) {
                if (coa[idx] > x) {
                    break;
                }
            }
            double left = idx == 0 ? 0 : coa[idx - 1];
            double right = idx < coa.length ? coa[idx] : 1;
            double r = (x - left) / (right - left);

            double thetaLeft = idx == 0 ? 0 : theta[idx - 1];
            double thetaRight = idx < theta.length ? theta[idx] : Math.PI / 2;

            return thetaLeft + r * (thetaRight - thetaLeft);

        }

        static void generateTable(String[] args) {

            for (double theta = Math.PI / 10; theta <= Math.PI / 2; theta += Math.PI / 10) {
                double num = 1 - Math.cos(theta);
                double denom = Math.sin(theta);
                double coa = num / denom;
                System.out.println(coa + "\t" + theta);
            }

        }
    }


}
