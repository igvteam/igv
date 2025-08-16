package org.broad.igv.renderer;

import org.broad.igv.track.RenderContext;

import java.awt.*;
import java.util.HashMap;
import java.util.Map;

public class DynSeqRenderer extends XYPlotRenderer {

    static Map<Character, Map<String, String>> letterPaths = new HashMap<>();
    Map<Character, Color> nucleotideColors = new HashMap<>();


    public DynSeqRenderer() {
        this.nucleotideColors = SequenceRenderer.getNucleotideColors();
    }

    @Override
    public String getDisplayName() {
        return "DynSeq";
    }

    /**
     * Render the data track as a bar chart.
     */
    @Override
    protected void drawDataPoint(Color graphColor, int dx, int pX, int baseY, int pY, RenderContext context) {

        int pixelsPerBP = (int) (1.0 / context.getReferenceFrame().getScale());

        if (pixelsPerBP < 5) {
            context.getGraphic2DForColor(graphColor).drawLine(pX, baseY, pX, pY);
        } else {
            renderDynSeq(context.getGraphics(), pixelsPerBP, pX, baseY, pY, 'C');
        }
        //}
    }


    private void renderDynSeq(Graphics2D ctx, int dx, int pX, int baseY, int pY, char base) {

        // Calculate rectangle position and height based on the wig value
        int rectY = Math.min(baseY, pY);
        int rectHeight = Math.abs(pY - baseY);

        boolean isNegative = pY > baseY;

        // Get nucleotide color from browser's color scheme
        Color nucleotideColor = nucleotideColors.getOrDefault(base, Color.GRAY);

        // Draw the base as a letter-shaped glyph
        drawLetterGlyph(ctx, base, pX, rectY, (int) dx, rectHeight, nucleotideColor, isNegative);


        // Draw overflow indicators if needed
//        if (feature.getValue() > dataRange.getMax()) {
//            IGVGraphics.fillRect(ctx, x, 0, width, 3, overflowColor);
//        } else if (feature.getValue() < dataRange.getMin()) {
//            IGVGraphics.fillRect(ctx, x, pixelHeight - 2, width, 3, overflowColor);
//        }
    }


    private void drawLetterGlyph(Graphics2D ctx, char base, int x, int y, int width, int height, Color color, boolean flipVertical) {
        Map<String, String> pathData = letterPaths.getOrDefault(base, letterPaths.get('N'));

        Graphics2D g2 = (Graphics2D) ctx.create();
        g2.setColor(color);

        if (flipVertical) {
            g2.translate(x + width / 2, y + height / 2);
            g2.scale(1, -1);
            g2.translate(-(x + width / 2), -(y + height / 2));
        }

        // Draw main path
        drawSVGPath(g2, pathData.get("main"), x, y, width, height);

        // Draw overlay path (if exists) - typically a cutout or highlight
        if (pathData.containsKey("overlay")) {
            g2.setColor(Color.WHITE);
            drawSVGPath(g2, pathData.get("overlay"), x, y, width, height);
        }

        g2.dispose();
    }

    private void drawSVGPath(Graphics2D ctx, String pathString, int x, int y, int width, int height) {
        // Parse SVG path string and draw it scaled to the given dimensions
        // Path is defined in 100x100 coordinate system
        double scaleX = width / 100.0;
        double scaleY = height / 100.0;

        // Enhanced SVG path parser for M (move), L (line), and C (cubic Bézier curve) commands
        java.util.regex.Pattern pattern = java.util.regex.Pattern.compile("[MLC][^MLC]*");
        java.util.regex.Matcher matcher = pattern.matcher(pathString);

        java.awt.geom.Path2D.Double path = new java.awt.geom.Path2D.Double();

        while (matcher.find()) {
            String command = matcher.group();
            char type = command.charAt(0);
            String[] coordStrs = command.substring(1).trim().split("[\\s,]+");
            double[] coords = new double[coordStrs.length];
            for (int i = 0; i < coordStrs.length; i++) {
                coords[i] = Double.parseDouble(coordStrs[i]);
            }

            if (type == 'M' && coords.length >= 2) {
                path.moveTo(x + coords[0] * scaleX, y + coords[1] * scaleY);
            } else if (type == 'L' && coords.length >= 2) {
                path.lineTo(x + coords[0] * scaleX, y + coords[1] * scaleY);
            } else if (type == 'C' && coords.length >= 6) {
                path.curveTo(
                        x + coords[0] * scaleX, y + coords[1] * scaleY, // control point 1
                        x + coords[2] * scaleX, y + coords[3] * scaleY, // control point 2
                        x + coords[4] * scaleX, y + coords[5] * scaleY  // end point
                );
            }
        }

        ctx.draw(path);
    }


    static {
        Map<String, String> aPaths = new HashMap<>();
        aPaths.put("main", "M 0 100 L 33 0 L 66 0 L 100 100 L 75 100 L 66 75 L 33 75 L 25 100 L 0 100");
        aPaths.put("overlay", "M 41 55 L 50 25 L 58 55 L 41 55");
        letterPaths.put('A', aPaths);

        Map<String, String> cPaths = new HashMap<>();
        cPaths.put("main", "M 100 28 C 100 -13 0 -13 0 50 C 0 113 100 113 100 72 L 75 72 C 75 90 30 90 30 50 C 30 10 75 10 75 28 L 100 28");
        letterPaths.put('C', cPaths);

        Map<String, String> gPaths = new HashMap<>();
        gPaths.put("main", "M 100 28 C 100 -13 0 -13 0 50 C 0 113 100 113 100 72 L 100 48 L 55 48 L 55 72 L 75 72 C 75 90 30 90 30 50 C 30 10 75 5 75 28 L 100 28");
        letterPaths.put('G', gPaths);

        Map<String, String> tPaths = new HashMap<>();
        tPaths.put("main", "M 0 0 L 0 20 L 35 20 L 35 100 L 65 100 L 65 20 L 100 20 L 100 0 L 0 0");
        letterPaths.put('T', tPaths);

        Map<String, String> nPaths = new HashMap<>();
        nPaths.put("main", "M 0 100 L 0 0 L 20 0 L 80 75 L 80 0 L 100 0 L 100 100 L 80 100 L 20 25 L 20 100 L 0 100");
        letterPaths.put('N', nPaths);
    }

}