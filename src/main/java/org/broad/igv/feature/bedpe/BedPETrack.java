package org.broad.igv.feature.bedpe;

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
import java.awt.geom.Arc2D;
import java.util.*;
import java.util.List;

import static org.broad.igv.feature.bedpe.BedPETrack.Direction.DOWN;
import static org.broad.igv.feature.bedpe.BedPETrack.Direction.UP;


/**
 * Created by jrobinso on 6/29/18.
 */
public class BedPETrack extends AbstractTrack {

    enum Direction {UP, DOWN}
    enum GraphType {block, arc}

    private  Genome genome;
    BedPETrack.Direction direction = UP; //DOWN;
    GraphType graphType =GraphType.arc;  // GraphType.block; //
    int thickness = 1;
    double logMaxScore = 0;

    private Map<String, List<BedPEFeature>> featureMap;
    private PEArcRenderer arcRenderer;
    private PEBLockRenderer blockRenderer;

    private int gap = 5;


    public BedPETrack() {
    }

    public BedPETrack(ResourceLocator locator, List<BedPEFeature> featureList, Genome genome) {
        super(locator);
        init(featureList, genome);
        this.genome = genome;
        arcRenderer = new PEArcRenderer();
        blockRenderer = new PEBLockRenderer();
        setHeight(250, true);
        setColor(new Color(180, 25, 137));
    }

    private void init(List<BedPEFeature> featureList, Genome genome) {

        this.featureMap = new HashMap<>();

        double maxScore = 0;
        for (BedPEFeature f : featureList) {
            String key = genome == null ? f.chr1 : genome.getCanonicalChrName(f.chr1);
            addToMap(f, key);
            if (!f.chr1.equals(f.chr2)) {
                key = genome == null ? f.chr2 : genome.getCanonicalChrName(f.chr2);
                addToMap(f, key);
            }

            maxScore = Math.max(maxScore, f.score);
        }

        if(maxScore > 0) logMaxScore = Math.log10(maxScore);


        featureMap.put(Globals.CHR_ALL, createWGFeatures(featureList, genome));

        // Sort feature lists by "start" (minimum of start1, start2)
        featureMap.values().forEach(flist -> Collections.sort(flist, (o1, o2) -> o1.getStart() - o2.getStart()));
    }

    private void addToMap(BedPEFeature f, String key) {
        List<BedPEFeature> features = featureMap.get(key);
        if (features == null) {
            features = new ArrayList<>();
            featureMap.put(key, features);
        }
        features.add(f);
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
            wgFeature.start2 = genome.getGenomeCoordinate(f.chr2, f.start2);
            wgFeature.end2 = genome.getGenomeCoordinate(f.chr2, f.end2);
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
    public void render(RenderContext context, Rectangle trackRectangle) {

        Graphics2D g2d = context.getGraphics();
        Rectangle clip = new Rectangle(g2d.getClip().getBounds());
        g2d.setClip(trackRectangle.intersection(clip.getBounds()));
        context.clearGraphicsCache();

        try {
            String chr = context.getReferenceFrame().getChrName();
            List<BedPEFeature> features = featureMap.get(chr);

            if (features != null) {

                Graphics2D g = (Graphics2D) context.getGraphics().create();
                double origin = context.getOrigin();
                double locScale = context.getScale();

                try {
                    for (BedPEFeature feature : features) {

                        double pixelStart = ((feature.getStart() - origin) / locScale);
                        double pixelEnd = ((feature.getEnd() - origin) / locScale);

                        // If the any part of the feature fits in the Track rectangle draw it
                        if (pixelEnd >= trackRectangle.getX() && pixelStart <= trackRectangle.getMaxX()) {
                            if (graphType == GraphType.block) {
                                blockRenderer.render(feature, context, trackRectangle, g);
                            } else {
                                arcRenderer.render(feature, context, trackRectangle, g);
                            }
                        }
                    }
                } finally {
                    g.dispose();
                }


            }
            context.clearGraphicsCache();
        } finally {
            g2d.setClip(clip);
        }
    }


    @Override
    public IGVPopupMenu getPopupMenu(TrackClickEvent te) {

        IGVPopupMenu menu = new IGVPopupMenu();

        menu.add(TrackMenuUtils.getTrackRenameItem(Collections.singleton(BedPETrack.this)));

        JMenuItem item = new JMenuItem("Set Track Height...");
        item.addActionListener(evt -> TrackMenuUtils.changeTrackHeight(Collections.singleton(BedPETrack.this)));
        menu.add(item);

        item = new JMenuItem("Set Track Color...");
        item.addActionListener(evt -> TrackMenuUtils.changeTrackColor(Collections.singleton(BedPETrack.this)));
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


        public void render(BedPEFeature feature, RenderContext context, Rectangle trackRectangle, Graphics2D g) {

            // TODO check chromosomes

            double origin = context.getOrigin();
            double locScale = context.getScale();

            Color trackColor = BedPETrack.this.getColor();
            double pixelStart = ((feature.getStart() - origin) / locScale);
            double pixelEnd = ((feature.getEnd() - origin) / locScale);
            double width = pixelEnd - pixelStart;

            int s = (feature.start1 + feature.end1) / 2;
            int e = (feature.start2 + feature.end2) / 2;

            double ps = (s - origin) / locScale;
            double pe = (e - origin) / locScale;

            Color fcolor = feature.color == null ? trackColor : feature.color;
            //if (fcolor != null && width > trackRectangle.width) {
            //    fcolor = getAlphaColor(fcolor, 0.1f);
            //}
            if (fcolor != null) {
                g.setColor(fcolor);
            }
            if (feature.thickness > 1) {
                g.setStroke(new BasicStroke(feature.thickness));
            }

            drawArc2(g, trackRectangle, ps, pe, feature.score, fcolor);
        }


        private Color getAlphaColor(Color fcolor, float alpha) {
            Color ac = alphaColors.get(fcolor);
            if (ac == null) {
                float[] rgb = new float[3];
                rgb = fcolor.getRGBColorComponents(rgb);
                ac = new Color(rgb[0], rgb[1], rgb[2], alpha);
                alphaColors.put(fcolor, ac);
            }
            return ac;
        }

        private void drawArc(Graphics2D g, Rectangle trackRectangle, double x1, double x2) {

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

        private void drawArc2(Graphics2D g, Rectangle trackRectangle, double x1, double x2, double score, Color color) {

            double pixelStart = Math.min(x1, x2);
            double pixelEnd = Math.max(x1, x2);

            int h = trackRectangle.height - gap;
            double logMax = BedPETrack.this.logMaxScore;
            if(logMax > 0 && score > 0){
                h = (int) ((Math.log10(score) / logMax) * h);
            }

            if (thickness > 1) {
                g.setStroke(new BasicStroke(thickness));
            }

            int w = (int) (pixelEnd - pixelStart);
            if (w < 3) {
                w = 3;
                pixelStart--;
            }

            double y = direction == UP ?   gap + trackRectangle.y + trackRectangle.height - h : gap + trackRectangle.y - h;
            int angleSt = direction == UP ? 0 : 180;
            Arc2D.Double arcPath = new Arc2D.Double(
                    pixelStart,
                    y,
                    w,
                    2*h,
                    angleSt,
                    180,
                    Arc2D.OPEN
            );

            g.draw(arcPath);

            Color shadedColor = getAlphaColor(color, 0.05f);
            g.setColor(shadedColor);
            g.fill(arcPath);
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


    class PEBLockRenderer {


        public void render(BedPEFeature feature, RenderContext context, Rectangle trackRectangle, Graphics2D g) {

            String chr = context.getChr();
            Genome genome = BedPETrack.this.genome;
            double origin = context.getOrigin();
            double locScale = context.getScale();

            Color trackColor = BedPETrack.this.getColor();
            Color fcolor = feature.color == null ? trackColor : feature.color;
            if (fcolor != null) {
                g.setColor(fcolor);
            }

            final int h = 10;
            final int blockY = trackRectangle.y + trackRectangle.height - h;

            int ps1 = (int) ((feature.start1 - origin) / locScale);
            int pe1 = (int) ((feature.end1 - origin) / locScale);
            String chr1 = genome == null ? feature.chr1 : genome.getCanonicalChrName(feature.chr1);
            if(chr1.equals(chr)) {
                // Trim width if possible to insure a gap between blocks
                int w1 = Math.max(1, pe1 - ps1);
                if (w1 > 3) w1--;
                if (w1 > 5) ps1++;
                if (pe1 >= trackRectangle.getX() && ps1 <= trackRectangle.getMaxX()) {
                    g.fillRect(ps1, blockY, w1, 10);
                }
            }


            int ps2 = (int) ((feature.start2 - origin) / locScale);
            int pe2 = (int) ((feature.end2 - origin) / locScale);
            String chr2 = genome == null ? feature.chr2 : genome.getCanonicalChrName(feature.chr2);
            if(chr2.equals(chr)) {
                // Trim width if possible to insure a gap between blocks
                int w2 = Math.max(1, pe2 - ps2);
                if (w2 > 3) w2--;
                if (w2 > 5) ps2++;

                if (pe2 >= trackRectangle.getX() && ps2 <= trackRectangle.getMaxX()) {
                    g.fillRect(ps2, blockY, w2, 10);
                }
            }

            // connecting line
            if (feature.isSameChr()) {
                int pl1 = Math.min(pe1, pe2);
                int pl2 = Math.max(ps1, ps2);
                final int connectorY = blockY + h/2;
                g.drawLine(pl1, connectorY, pl2, connectorY);
            }
        }

    }

    @Override
    public void marshalXML(Document document, Element element) {

        super.marshalXML(document, element);

        element.setAttribute("direction", String.valueOf(direction));
        element.setAttribute("thickness", String.valueOf(thickness));
        element.setAttribute("graphType", String.valueOf(graphType));

    }

    @Override
    public void unmarshalXML(Element element, Integer version) {

        super.unmarshalXML(element, version);

        if (element.hasAttribute("direction"))
            this.direction = Direction.valueOf(element.getAttribute("direction"));
        if (element.hasAttribute("thickness"))
            this.thickness = Integer.parseInt(element.getAttribute("thickness"));
        if (element.hasAttribute("graphType"))
            this.graphType = GraphType.valueOf(element.getAttribute("graphType"));

    }
}
