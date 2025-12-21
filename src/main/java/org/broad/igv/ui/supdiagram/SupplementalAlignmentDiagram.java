package org.broad.igv.ui.supdiagram;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.sam.SupplementaryAlignment;
import org.broad.igv.sam.SupplementaryGroup;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.util.IGVMouseInputAdapter;
import org.broad.igv.util.ChromosomeColors;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.geom.CubicCurve2D;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

class SupplementalAlignmentDiagram extends JPanel {

    public static final int LABEL_GAP = 2;
    private static final Logger log = LogManager.getLogger(SupplementalAlignmentDiagram.class);

    public static final int BORDER_GAP = 30;
    public static final int BETWEEN_ALIGNMENT_GAP = 15;
    public static final int BETWEEN_CONTIG_GAP = 10;
    public static final int ALIGNMENT_HEIGHT = 10;
    public static final Color SELECTED_COLOR = Color.BLUE;
    public static final Color PRIMARY_BORDER_COLOR = Color.DARK_GRAY;
    public static final int DEFAULT_WIDTH = 500;

    private Rectangle2D chrDiagramBounds = null;
    private Rectangle2D readDiagramBounds = null;

    //This keeps track of the polygons on screen so it's possible to do bounds checks for mouse interaction
    private final Map<AlignmentArrow, SupplementaryAlignment> elementsOnScreen = new LinkedHashMap<>();
    private final Set<SupplementaryAlignment> selected = new LinkedHashSet<>();
    private final SupplementaryGroup toDraw;

    public SupplementalAlignmentDiagram(final SupplementaryGroup toDraw) {
        this.toDraw = toDraw;
        this.setBackground(Color.WHITE);
        this.addMouseMotionListener(new IGVMouseInputAdapter() {

            @Override
            public void mouseMoved(final MouseEvent e) {
                selected.clear();
                //dumb brute force search but there should only ever be a handful of these so it should be fine
                for (AlignmentArrow arrow : elementsOnScreen.keySet()) {
                    if (arrow.contains(e.getPoint())) {
                        selected.add(elementsOnScreen.get(arrow));
                        repaint();
                        return;
                    }
                }
                repaint();
            }
        });

        this.addMouseListener(new IGVMouseInputAdapter() {
            @Override
            public void mouseClicked(final MouseEvent e) {
                if (SwingUtilities.isLeftMouseButton(e)) {
                    IGV igv = null; // for testing we may not have an IGV instance available
                    try {
                        igv = IGV.getInstance();
                    } catch (RuntimeException ex) {
                        log.info("Clicked " + selected.stream().map(SupplementaryAlignment::toString)
                                .collect(Collectors.joining("\n")));
                    }

                    if (igv != null && !selected.isEmpty()) {
                        final SupplementaryAlignment first = selected.stream().findFirst().get();
                        igv.getAlignmentTracks().forEach(
                                t -> t.setSelectedAlignment(toDraw.unwrap())
                        );
                        if (e.isShiftDown()) { //if shift is held add to the existing set of frames instead of replacing them
                            FrameManager.addNewLociToFrames(FrameManager.getDefaultFrame(), java.util.List.of(first));
                        } else {
                            igv.setDefaultFrame(first.getContig() + ":" + first.getStart() + "-" + first.getEnd());
                        }
                    }
                }
            }
        });
    }


    @Override
    protected void paintComponent(final Graphics g) {
        super.paintComponent(g);
        setBackground(Color.WHITE);
        ((Graphics2D) g).setComposite(getAlphaComposite());
        g.setColor(Color.LIGHT_GRAY);
        elementsOnScreen.clear();
        final int halfHeight = getHeight() / 2;
        final int width = getWidth();
        chrDiagramBounds = new Rectangle2D.Float(0, 10, width, halfHeight - 10);
        drawContigOrder(g, chrDiagramBounds);
        readDiagramBounds = new Rectangle2D.Float(0, halfHeight, width, halfHeight);
        drawReadOrder(g, readDiagramBounds);

        if (!selected.isEmpty()) {
            final SupplementaryAlignment first = this.selected.iterator().next();
            g.setColor(Color.BLACK);
            g.drawString(String.format("%s:%d-%d", first.getContig(), first.getStart(), first.getEnd()), 30, getHeight() - g.getFontMetrics().getHeight() + 2);
        }
    }

    private void drawReadOrder(final Graphics g, final Rectangle2D bounds) {
        g.setColor(Color.DARK_GRAY);
        ((Graphics2D) g).draw(bounds);
        g.drawString("Read Order", (int) bounds.getX() + 5, (int) bounds.getY() + 15);
        Map<SupplementaryAlignment, AlignmentArrow> saInReadOrder = drawAlignmentsInReadOrder((Graphics2D) g.create(), toDraw, selected, bounds);
        drawArcs((Graphics2D) g.create(), toDraw, selected, saInReadOrder, bounds);
        final int lowestArrowPoint = saInReadOrder.values().stream().mapToInt(a -> (int) a.getBounds2D().getMaxY()).max().getAsInt();
        drawReadLengthLabel(((Graphics2D) g.create()), lowestArrowPoint + 15, saInReadOrder);
        saInReadOrder.forEach((k, v) -> elementsOnScreen.put(v, k));
        g.setColor(Color.DARK_GRAY);
        ((Graphics2D) g).draw(bounds);
    }

    private void drawContigOrder(final Graphics g, final Rectangle2D bounds) {
        g.setColor(Color.DARK_GRAY);
        ((Graphics2D) g).draw(bounds);
        g.drawString("Alignment Order", (int) bounds.getX() + 5, (int) bounds.getY() + 15);
        Map<SupplementaryAlignment, AlignmentArrow> saInPositionOrder = drawAlignmentsInCondensedChromosomeOrder((Graphics2D) g.create(), toDraw, selected, bounds);
        drawArcs((Graphics2D) g.create(), toDraw, selected, saInPositionOrder, bounds);
        final int lowestArrowPoint = saInPositionOrder.values().stream().mapToInt(a -> (int) a.getBounds2D().getMaxY()).max().getAsInt();
        drawContigLabels((Graphics2D) g.create(), lowestArrowPoint + 15, saInPositionOrder);
        saInPositionOrder.forEach((k, v) -> elementsOnScreen.put(v, k));
        g.setColor(Color.DARK_GRAY);
        ((Graphics2D) g).draw(bounds);
    }

    private void drawReadLengthLabel(final Graphics2D g, final int mid, final Map<SupplementaryAlignment, AlignmentArrow> saInReadOrder) {
        java.util.List<AlignmentArrow> arrowsInOrder = new ArrayList<>(saInReadOrder.values());
        int left = (int) arrowsInOrder.get(0).getBounds().getMinX();
        int right = (int) arrowsInOrder.get(arrowsInOrder.size() - 1).getBounds().getMaxX();
        final int baseCount = toDraw.getBaseCount();
        g.setColor(Color.BLACK);
        drawCenteredStringWithRangeLines(g, mid, "Length in bases = " + baseCount, left, right);
    }

    private static void drawArcs(final Graphics2D g,
                                 final SupplementaryGroup toDraw,
                                 final Set<SupplementaryAlignment> selected,
                                 final Map<SupplementaryAlignment, AlignmentArrow> saToArrowMap,
                                 final Rectangle2D bounds) {
        toDraw.iterateInReadOrder()
                .forEachRemaining(sa -> {
                    SupplementaryAlignment next = toDraw.getNextInRead(sa);
                    if (next != null) {
                        final boolean highlight = selected.contains(sa) || selected.contains(next);
                        drawArc(g, saToArrowMap.get(sa), saToArrowMap.get(next), highlight, bounds);
                    }
                });
    }

    private static Composite getAlphaComposite() {
        final float alpha = 0.75f;
        final int type = AlphaComposite.SRC_OVER;
        return AlphaComposite.getInstance(type, alpha);
    }

    private static Map<SupplementaryAlignment, AlignmentArrow> drawAlignmentsInReadOrder(final Graphics2D g, final SupplementaryGroup toDraw,
                                                                                         final Set<SupplementaryAlignment> selected, Rectangle2D bounds) {
        int midline = (int) (bounds.getY() + .5 * bounds.getHeight());
        final Map<SupplementaryAlignment, AlignmentArrow> positions = new LinkedHashMap<>();
        final int totalAlignedBases = toDraw.getBaseCount();
        final int scaledAlignmentGap = scale(2, BETWEEN_ALIGNMENT_GAP, bounds.getWidth());
        final int scaledBorderGap = scale(BORDER_GAP / 3, BORDER_GAP, bounds.getWidth());
        final double availableSpace = bounds.getWidth() - (2 * scaledBorderGap + (toDraw.size() - 1) * scaledAlignmentGap);

        double lastPosition = scaledBorderGap;
        for (SupplementaryAlignment sa : (Iterable<SupplementaryAlignment>) toDraw::iterateInReadOrder) {
            final double spaceToUse = getSpaceToUse(availableSpace, sa.getNumberOfAlignedBases(), totalAlignedBases);
            final int end = (int) (lastPosition + spaceToUse);
            AlignmentArrow readArrow = new AlignmentArrow(midline, ALIGNMENT_HEIGHT, (int) (bounds.getX() + lastPosition), (int) (bounds.getX() + end), sa.getStrand());
            lastPosition = end + scaledAlignmentGap;
            positions.put(sa, readArrow);
        }

        drawArrows(g, selected, positions, toDraw.getPrimaryAlignment());
        return positions;
    }

    private static Graphics2D getSelectedGraphics(final Graphics2D g) {
        Graphics2D selectedGraphics = (Graphics2D) g.create();
        selectedGraphics.setStroke(new BasicStroke(3.0f));
        selectedGraphics.setColor(SELECTED_COLOR);
        return selectedGraphics;
    }

    private static <T extends Locatable> Map<Locatable, List<T>> groupBySpanningInterval(List<T> intervals) {
        List<T> currentGroup = null;
        Map<Locatable, List<T>> output = new LinkedHashMap<>();
        Locatable spanning = null;
        if (intervals.isEmpty()) {
            return Collections.emptyMap();
        }
        for (T loc : intervals) {
            if (currentGroup == null || currentGroup.isEmpty()) {
                currentGroup = new ArrayList<>();
                currentGroup.add(loc);
                spanning = loc;
            } else if (spanning.overlaps(loc)) {
                currentGroup.add(loc);
                spanning = new Interval(spanning.getContig(), Math.min(spanning.getStart(), loc.getStart()), Math.max(spanning.getEnd(), loc.getEnd()));
            } else {
                output.put(spanning, currentGroup);
                currentGroup = new ArrayList<>();
                currentGroup.add(loc);
                spanning = loc;
            }
        }
        output.put(spanning, currentGroup);
        return output;
    }

    //Draw the alignments in chromosome order but give an indication of how close they are to each other i.e. within a 1kb window or not
    //Handle overlapping alignments
    private static Map<SupplementaryAlignment, AlignmentArrow> drawAlignmentsInCondensedChromosomeOrder(final Graphics2D g, final SupplementaryGroup toDraw,
                                                                                                        final Set<SupplementaryAlignment> selected, Rectangle2D bounds) {


        final double midline = bounds.getY() + .5 * bounds.getHeight();
        final Map<SupplementaryAlignment, AlignmentArrow> positions = new LinkedHashMap<>();
        final List<String> contigs = toDraw.getContigs();
        final int scaledAlignmentGap = scale(2, BETWEEN_ALIGNMENT_GAP, bounds.getWidth());
        final int scaledContigGap = scale(2, BETWEEN_CONTIG_GAP, bounds.getWidth());
        final int scaledBorderGap = scale(BORDER_GAP / 3, BORDER_GAP, bounds.getWidth());

        final var groupedBySpan = groupBySpanningInterval(toDraw.streamInPositionOrder().toList());
        final Map<String, List<Locatable>> byContig = groupedBySpan.keySet()
                .stream()
                .collect(Collectors.groupingBy(Locatable::getContig, LinkedHashMap::new, Collectors.toList()));


        //tihs should probably vary per contig instead of being uniform
        final double perContigAvailableSpace = (bounds.getWidth() - (2 * scaledBorderGap + (contigs.size() - 1) * scaledContigGap)) / ((double) contigs.size());

        int contigStart = scaledBorderGap;
        for (var contigEntry : byContig.entrySet()) {
            //find the available space for each contig and set the drawing head there
            int contigEnd = (int) (contigStart + perContigAvailableSpace);
            List<Locatable> distinctSpans = contigEntry.getValue();
            //find the reference length of all the span groups on this contig
            int totalSpansLength = distinctSpans.stream().mapToInt(Locatable::getLengthOnReference).sum();
            int spanStart = contigStart;
            for (Locatable span : distinctSpans) {
                //now handle each span group
                int spanLength = span.getLengthOnReference();
                int spanSpaceAvailable = (int) getSpaceToUse((double) perContigAvailableSpace - (distinctSpans.size() - 1) * scaledAlignmentGap, spanLength, totalSpansLength);
                final List<SupplementaryAlignment> activeSpanGroup = groupedBySpan.get(span);
                for (int i = 0; i < activeSpanGroup.size(); i++) {
                    //each read in the span is placed relatively within the space
                    SupplementaryAlignment sa = activeSpanGroup.get(i);
                    int scaledReadStart = (int) getSpaceToUse(spanSpaceAvailable, sa.getStart() - span.getStart(), spanLength);
                    int scaledReadEnd = (int) getSpaceToUse(spanSpaceAvailable, sa.getEnd() - span.getStart(), spanLength);
                    // height offset is used to try to layout overlapping arrows above instead of ontop of each other
                    int heightOffset = ALIGNMENT_HEIGHT - (int) ((2 * ALIGNMENT_HEIGHT) * (1.0 / (activeSpanGroup.size() + 1.0)) * (i + 1.0));
                    final AlignmentArrow readArrow = new AlignmentArrow((int) midline + 2 * heightOffset, ALIGNMENT_HEIGHT,
                            (int) bounds.getX() + spanStart + scaledReadStart,
                            (int) bounds.getX() + spanStart + scaledReadEnd,
                            sa.getStrand());

                    positions.put(sa, readArrow);

                }
                spanStart += spanSpaceAvailable + scaledAlignmentGap;
            }
            //move contig start forward
            contigStart = contigEnd + scaledContigGap;
        }

        // | ... [  ]>-5kbp<[  ] ... [  ]> |     |[   |>  ... |
        // Contigs all the same scale?
        // Contigs proportional to actual size?  <
        //1 divide available space into contig zones
        //2 go through all reads on contig and count necessary zones.
        //  edges if the read hits left / right edge
        // merge overlappers into single zones
        // discover close by / far away zones

        drawArrows(g, selected, positions, toDraw.getPrimaryAlignment());
        return positions;
    }

    //This draws the aligments scaled according to their aligned length on the reference and in chromosome order
    //TODO fix overlapping alignments
    private static Map<SupplementaryAlignment, AlignmentArrow> drawAlignmentsInPositionOrder(final Graphics2D g, final SupplementaryGroup toDraw,
                                                                                             final Set<SupplementaryAlignment> selected, final int width, final int midline) {
        final Map<SupplementaryAlignment, AlignmentArrow> positions = new LinkedHashMap<>();
        final List<String> contigs = toDraw.getContigs();
        final int totalAlignedBases = toDraw.getLengthOnReference();
        final int scaledAlignmentGap = scale(2, BETWEEN_ALIGNMENT_GAP, width);
        final int scaledContigGap = scale(2, BETWEEN_CONTIG_GAP, width);
        final int scaledBorderGap = scale(BORDER_GAP / 3, BORDER_GAP, width);
        final double availableSpace = width - (2 * scaledBorderGap + (toDraw.size() - 1) * scaledAlignmentGap + (contigs.size() - 1) * scaledContigGap);

        String lastContig = contigs.get(0);
        double lastPosition = scaledBorderGap;
        for (SupplementaryAlignment sa : (Iterable<SupplementaryAlignment>) toDraw::iterateInPositionOrder) {

            final double spaceToUse = getSpaceToUse(availableSpace, sa.getLengthOnReference(), totalAlignedBases);
            final String newContig = sa.getContig();
            if (lastPosition != scaledBorderGap && !Objects.equals(lastContig, newContig)) {
                lastPosition += scaledContigGap;
            }

            lastContig = newContig;
            final int end = (int) (lastPosition + spaceToUse);

            AlignmentArrow readArrow = new AlignmentArrow(midline, ALIGNMENT_HEIGHT, (int) lastPosition, end, sa.getStrand());
            positions.put(sa, readArrow);
            lastPosition = end + scaledAlignmentGap;
        }

        drawArrows(g, selected, positions, toDraw.getPrimaryAlignment());
        return positions;
    }

    /**
     * Scale a value between min/max based on the difference between width and DEFAULT_WIDTH
     */
    private static int scale(final int min, final int max, final double width) {
        if (width >= DEFAULT_WIDTH) {
            return max;
        } else {
            final double scaleDown = width / (double) DEFAULT_WIDTH;
            return Math.max((int) (max * scaleDown), min);
        }
    }

    private static void drawArrows(final Graphics2D g, final Set<SupplementaryAlignment> selected, final Map<SupplementaryAlignment, AlignmentArrow> positions, final SupplementaryAlignment primary) {
        //first draw background colors
        positions.forEach((sa, readArrow) -> {
            g.setColor(ChromosomeColors.getColor(sa.getContig()));
            g.fill(readArrow);
        });

        //outline primary read
        g.setColor(PRIMARY_BORDER_COLOR);
        final AlignmentArrow primaryShape = positions.get(primary);
        g.draw(primaryShape);

        //next draw borders on top to work in case of overlaps, this should come last so it doesn't get clobbered
        //in the case of overlaps
        for (Map.Entry<SupplementaryAlignment, AlignmentArrow> pair : positions.entrySet()) {
            SupplementaryAlignment sa = pair.getKey();
            AlignmentArrow readArrow = pair.getValue();
            g.setColor(ChromosomeColors.getColor(sa.getContig()));
            Graphics2D g2 = selected.contains(sa) ? getSelectedGraphics(g) : g;
            g2.draw(readArrow);
        }
    }

    private static double getSpaceToUse(final double availableSpace, final int numberOfBasePairs, final int totalAlignedBases) {
        final double fractionOfWhole = (double) numberOfBasePairs / totalAlignedBases;
        return fractionOfWhole * availableSpace;
    }

    private static void drawContigLabels(final Graphics2D g, final int mid, final Map<SupplementaryAlignment, AlignmentArrow> positions) {
        positions.keySet().
                stream()
                .collect(Collectors.groupingBy(SupplementaryAlignment::getContig))
                .forEach((contig, alignments) -> {
                    int start = (int) positions.get(alignments.get(0)).getBounds().getX();
                    final Rectangle rightmostBounds = positions.get(alignments.get(alignments.size() - 1)).getBounds();
                    int end = (int) (rightmostBounds.getX() + rightmostBounds.getWidth());
                    g.setColor(ChromosomeColors.getColor(contig));
                    drawCenteredStringWithRangeLines(g, mid, contig, start, end);
                });
    }

    private static void drawCenteredStringWithRangeLines(final Graphics2D g, final int mid, final String string, final int start, final int end) {
        final FontMetrics fontMetrics = g.getFontMetrics();
        final int labelWidth = fontMetrics.stringWidth(string);
        //    s----label-----e
        final int lineY = mid;
        final int height = fontMetrics.getHeight();
        if (labelWidth + 2 * LABEL_GAP < end - start) {
            final int leftLineEnd = (start + end - labelWidth) / 2 - LABEL_GAP;
            final int rightLineStart = (start + end + labelWidth) / 2 + LABEL_GAP;
            g.drawLine(start, lineY - 2, start, lineY + 2);
            g.drawLine(end, lineY - 2, end, lineY + 2);
            g.drawLine(start, lineY, leftLineEnd, lineY);
            g.drawString(string, leftLineEnd + LABEL_GAP, lineY + height / 3);
            g.drawLine(rightLineStart, lineY, end, lineY);
        } else {
            g.drawString(string, start, lineY + height / 3);
        }
    }

    private static void drawArc(final Graphics2D g, final AlignmentArrow currentPos, final AlignmentArrow nextPos, final boolean highlight, Rectangle2D bounds) {
        CubicCurve2D c = new CubicCurve2D.Double();
        Point2D from = currentPos.getTip();
        Point2D to = nextPos.getTail();
        g.setStroke(new BasicStroke(2, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER));


        // Determine the distance between the two points
        double distance = Math.abs(from.getX() - to.getX());

        // Calculate the height proportional to the distance. We'll use 1/4 of the bounding box height
        // as a baseline and then adjust based on the distance.
        double minHeight = bounds.getY();
        double heightAdjustment = 1; //.75 + .25 * (distance / bounds.getWidth());
        double maxHeight = bounds.getHeight() - 5;
        double actualHeight = maxHeight * heightAdjustment;
        double controlY = bounds.getY() + (bounds.getHeight() - actualHeight);

        // Create two control points for the Bezier curve.
        Point2D control1 = new Point2D.Double(from.getX() + (to.getX() - from.getX()) / 4, controlY);
        Point2D control2 = new Point2D.Double(from.getX() + 3 * (to.getX() - from.getX()) / 4, controlY);


        c.setCurve(from, control1, control2, to);
        g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        if (highlight) {
            g.setColor(SupplementaryAlignmentDiagramDialog.ARC_HIGHLIGHT_COLOR);
        } else {
            g.setColor(Color.GRAY);
        }
        g.draw(c);
        drawPoint(g, to);
    }


    private static void drawPoint(Graphics2D g, Point2D point) {
        g.draw(new Ellipse2D.Double(point.getX() - 1, point.getY() - 1, 2, 2));
    }

    private static void writeContigName(final Graphics g, final String contig, final double x, final int y) {
        Color originalColor = g.getColor();
        try {
            g.setColor(ChromosomeColors.getColor(contig));
            g.drawString(contig, (int) x, y);
        } finally {
            g.setColor(originalColor);
        }
    }
}
