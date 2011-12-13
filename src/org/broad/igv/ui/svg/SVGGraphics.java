/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.ui.svg;

import org.broad.igv.ui.color.ColorUtilities;

import javax.swing.*;
import java.awt.*;
import java.awt.font.FontRenderContext;
import java.awt.font.GlyphVector;
import java.awt.geom.AffineTransform;
import java.awt.image.BufferedImage;
import java.awt.image.BufferedImageOp;
import java.awt.image.ImageObserver;
import java.awt.image.RenderedImage;
import java.awt.image.renderable.RenderableImage;
import java.io.PrintWriter;
import java.text.AttributedCharacterIterator;
import java.util.Map;

/**
 * This is an experimental class, its not actually used (yet)
 *
 * @date Nov 20, 2010
 */
public class SVGGraphics extends Graphics2D {

    PrintWriter outputStream;
    Color color = Color.black;
    Color background = Color.white;
    BasicStroke stroke;
    Font font;
    FontRenderContext fontRenderContext;
    FontMetrics fontMetrics;
    AffineTransform transform;
    JComponent component; // <= dummy component to obtain font metrics
    int translateX;
    int translateY;
    String clipRef;
    static int clipNumber = 0;
    Shape clip;


    public SVGGraphics(PrintWriter outputStream, Rectangle clip) {
                   this.clip = clip;
        initDefaults();

        this.outputStream = outputStream;
    }

    public SVGGraphics(SVGGraphics g) {
        this.clip = g.clip;
        this.outputStream = g.outputStream;
        this.color = g.color;
        this.background = g.background;
        this.stroke = g.stroke;
        this.font = g.font;
        this.fontMetrics = g.fontMetrics;
        this.fontRenderContext = g.fontRenderContext;
        this.component = g.component;
        if (g.transform != null) {
            this.transform = new AffineTransform(g.transform);
        }
        this.clipRef = g.clipRef;
        this.translateX = g.translateX;
        this.translateY = g.translateY;
    }


    private void initDefaults() {
        component = new JPanel();
        stroke = new BasicStroke();
        color = Color.black;
        font = new Font("Arial", Font.PLAIN, 12);
        fontMetrics = component.getFontMetrics(font);
        fontRenderContext = new FontRenderContext(null, true, false);
        translateX = 0;
        translateY = 0;
    }

    @Override
    public void draw(Shape shape) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public boolean drawImage(Image image, AffineTransform affineTransform, ImageObserver imageObserver) {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void drawImage(BufferedImage bufferedImage, BufferedImageOp bufferedImageOp, int i, int i1) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void drawRenderedImage(RenderedImage renderedImage, AffineTransform affineTransform) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void drawRenderableImage(RenderableImage renderableImage, AffineTransform affineTransform) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void drawString(String s, int x, int y) {

        int tx = translateX + x;
        int ty = translateY + y;

        outputStream.print("<text x=\"");
        outputStream.print(tx);
        outputStream.print("\" y=\"");
        outputStream.print(ty);
        outputStream.println("\"");
        outputStream.print(" style=\"font-family: " + font.getName() + ";");
        outputStream.print(" font-size:" + font.getSize() + ";");

        if (font.isBold()) {
            outputStream.print(" font-weight: bold;");
        }

        outputStream.print(" fill:rgb" + getRGBString());
        
        if (clipRef != null) {
            outputStream.print(";clip-path: " + clipRef);
        }

        outputStream.print("\" ");

        applyTransform();
        
        outputStream.print(">");
        outputStream.print(s);
        outputStream.println("</text>");
    }

    @Override
    public void drawString(String s, float v, float v1) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void drawString(AttributedCharacterIterator attributedCharacterIterator, int i, int i1) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public boolean drawImage(Image image, int i, int i1, ImageObserver imageObserver) {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public boolean drawImage(Image image, int i, int i1, int i2, int i3, ImageObserver imageObserver) {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public boolean drawImage(Image image, int i, int i1, Color color, ImageObserver imageObserver) {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public boolean drawImage(Image image, int i, int i1, int i2, int i3, Color color, ImageObserver imageObserver) {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public boolean drawImage(Image image, int i, int i1, int i2, int i3, int i4, int i5, int i6, int i7, ImageObserver imageObserver) {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public boolean drawImage(Image image, int i, int i1, int i2, int i3, int i4, int i5, int i6, int i7, Color color, ImageObserver imageObserver) {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void dispose() {
    }

    @Override
    public void drawString(AttributedCharacterIterator attributedCharacterIterator, float v, float v1) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void drawGlyphVector(GlyphVector glyphVector, float v, float v1) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void fill(Shape shape) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public boolean hit(Rectangle rectangle, Shape shape, boolean b) {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public GraphicsConfiguration getDeviceConfiguration() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void setComposite(Composite composite) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void setPaint(Paint paint) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void setStroke(Stroke stroke) {
        if (!(stroke instanceof BasicStroke)) {
            // throw exception
        }
        this.stroke = (BasicStroke) stroke;
    }

    @Override
    public void setRenderingHint(RenderingHints.Key key, Object o) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Object getRenderingHint(RenderingHints.Key key) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void setRenderingHints(Map<?, ?> map) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void addRenderingHints(Map<?, ?> map) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public RenderingHints getRenderingHints() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Graphics create() {
        return new SVGGraphics(this);
    }

    @Override
    public void translate(int x, int y) {
        translateX += x;
        translateY += y;
    }


    @Override
    public void translate(double x, double y) {
        System.out.println("Translate " + x + " " + y);
    }
    

    @Override
    public Color getColor() {
        return color;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void setColor(Color color) {
        this.color = color;
    }

    @Override
    public void setPaintMode() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void setXORMode(Color color) {
        //To change body of implemented methods use File | Settings | File Templates.
    }


    @Override
    public Font getFont() {
        return font;
    }

    @Override
    public void setFont(Font font) {
        this.font = font;
        fontMetrics = component.getFontMetrics(font);
    }

    @Override
    public FontMetrics getFontMetrics() {
        return fontMetrics;
    }

    @Override
    public FontMetrics getFontMetrics(Font font) {
        return component.getFontMetrics(font);
    }

    @Override
    public Rectangle getClipBounds() {
        return clip.getBounds();
    }

    @Override
    public void clipRect(int x, int y, int w, int h) {
        setClip(x, y, w, h);
    }


    @Override
    public void setClip(int x, int y, int w, int h) {

            clip = new Rectangle(x, y, w, h);
            int tx = translateX + x;
            int ty = translateY + y;
            String id = getClipId();
            outputStream.println("<clipPath id=\"" + id + "\">\n" +
                    "    <rect id=\"rect1\" x=\"" + tx + "\" y=\"" + ty + "\"\n" +
                    "       width=\"" + w + "\" height=\"" + h + "\"/>");
            outputStream.println("</clipPath>");

            clipRef = "url(#" + id + ")";
    }

    private static synchronized String  getClipId() {
          return "clip_" + clipNumber++;
    }




    @Override
    public Shape getClip() {
        return clip;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void setClip(Shape shape) {
        this.clip = shape;
    }

    @Override
    public void copyArea(int i, int i1, int i2, int i3, int i4, int i5) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void drawLine(int x1, int y1, int x2, int y2) {

        int tx1 = translateX + x1;
        int ty1 = translateY + y1;
        int tx2 = translateX + x2;
        int ty2 = translateY + y2;
        String rgbString = getRGBString();
        float strokeWidth = stroke.getLineWidth();

        outputStream.println("<line x1=\"" + tx1 + "\" y1=\"" + ty1 + "\" x2=\"" + tx2 + "\" y2=\"" + ty2 + "\"");
        outputStream.print("style=\"stroke:rgb" + rgbString + ";stroke-width:" + strokeWidth);

        if (clipRef != null) {
            outputStream.print(";clip-path: " + clipRef);
        }

        applyTransform();

        outputStream.println("\"/>");
    }

    @Override
    public void drawRect(int x, int y, int w, int h) {

        int tx = translateX + x;
        int ty = translateY + y;
        String rgbString = getRGBString();
        float strokeWidth = stroke.getLineWidth();

        outputStream.print("<rect x=\"");
        outputStream.print(tx);
        outputStream.print("\" y=\"");
        outputStream.print(ty);
        outputStream.print("\" width=\"");
        outputStream.print(w);
        outputStream.print("\" height=\"");
        outputStream.print(h);
        outputStream.println("\"");
        outputStream.print("  style=\"fill:none;stroke-width:" + strokeWidth + ";stroke:rgb" + rgbString);

        if (clipRef != null) {
            outputStream.print(";clip-path: " + clipRef);
        }

        applyTransform();

        outputStream.println("\"/>");
    }

    private String getRGBString() {
        String rgbString = ColorUtilities.colorToString(color);
        return "(" + rgbString + ")";
    }

    @Override
    public void fillRect(int x, int y, int w, int h) {

        int tx = translateX + x;
        int ty = translateY + y;
        String rgbString = getRGBString();
        float strokeWidth = stroke.getLineWidth();

        outputStream.print("<rect x=\"");
        outputStream.print(tx);
        outputStream.print("\" y=\"");
        outputStream.print(ty);
        outputStream.print("\" width=\"");
        outputStream.print(w);
        outputStream.print("\" height=\"");
        outputStream.print(h);
        outputStream.println("\"");
        outputStream.print("  style=\"fill:rgb" + rgbString);

        if (clipRef != null) {
            outputStream.print(";clip-path: " + clipRef);
        }

        applyTransform();

        outputStream.println("\"/>");

    }


    @Override
    public void drawPolygon(int[] x, int[] y, int npts) {
        String rgbString = getRGBString();
        float strokeWidth = stroke.getLineWidth();

        outputStream.print("<polygon points=\"");
        for (int i = 0; i < npts; i++) {
            int tx = translateX + x[i];
            int ty = translateY + y[i];
            outputStream.print(tx);
            outputStream.print(",");
            outputStream.print(ty);
            outputStream.print(" ");
        }
        outputStream.println("\"");
        outputStream.print("  style=\"fill:none;stroke-width:" + strokeWidth + ";stroke:rgb" + rgbString);

        if (clipRef != null) {
            outputStream.print(";clip-path: " + clipRef);
        }

        applyTransform();

        outputStream.println("\"/>");
    }

    @Override
    public void fillPolygon(int[] x, int[] y, int npts) {
        String rgbString = getRGBString();
        float strokeWidth = stroke.getLineWidth();

        outputStream.print("<polygon points=\"");
        for (int i = 0; i < npts; i++) {
            int tx = translateX + x[i];
            int ty = translateY + y[i];
            outputStream.print(tx);
            outputStream.print(",");
            outputStream.print(ty);
            outputStream.print(" ");
        }
        outputStream.println("\"");
        outputStream.print("  style=\"fill:rgb" + rgbString);

        if (clipRef != null) {
            outputStream.print(";clip-path: " + clipRef);
        }

        applyTransform();

        outputStream.println("\"/>");
    }

    @Override
    public void clearRect(int x, int y, int w, int h) {
        Color c = color;
        try {
            color = background;
            fillRect(x, y, w, h);
        }
        finally {
            color = c;
        }
    }

    @Override
    public void drawRoundRect(int i, int i1, int i2, int i3, int i4, int i5) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void fillRoundRect(int i, int i1, int i2, int i3, int i4, int i5) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void drawOval(int i, int i1, int i2, int i3) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void fillOval(int i, int i1, int i2, int i3) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void drawArc(int i, int i1, int i2, int i3, int i4, int i5) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void fillArc(int i, int i1, int i2, int i3, int i4, int i5) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void drawPolyline(int[] ints, int[] ints1, int i) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void rotate(double v) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void rotate(double v, double v1, double v2) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void scale(double v, double v1) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void shear(double v, double v1) {
        //To change body of implemented methods use File | Settings | File Templates.
    }


    @Override
    public void transform(AffineTransform affineTransform) {
        if (transform == null) {
            setTransform(affineTransform);
        } else {
            transform.concatenate(affineTransform);
        }
    }

    @Override
    public void setTransform(AffineTransform affineTransform) {
        this.transform = affineTransform;
    }


    private void applyTransform() {
        if (transform != null) {

            // TODO -- cache the matrix
            double[] matrix = new double[6];
            transform.getMatrix(matrix);
            outputStream.print(" transform=\"matrix(");
            for (int i = 0; i < matrix.length; i++) {
                outputStream.print(matrix[i]);
                if (i < matrix.length - 1) outputStream.print(" ");
            }
            outputStream.println(")\" ");


        }
    }

    @Override
    public AffineTransform getTransform() {
        return transform;
    }

    @Override
    public Paint getPaint() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Composite getComposite() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void setBackground(Color color) {
        this.background = color;
    }

    @Override
    public Color getBackground() {
        return background;  //To change body of implemented methods use File | Settings | File Templates.
    }


    @Override
    public Stroke getStroke() {
        return stroke;  //To change body of implemented methods use File | Settings | File Templates.
    }


    @Override
    public void clip(Shape shape) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public FontRenderContext getFontRenderContext() {
        return fontRenderContext;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
