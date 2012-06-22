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
/**
 * SnapshotUtilities.java
 *
 * Created on November 29, 2007, 2:14 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package org.broad.igv.ui.util;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.apache.log4j.Logger;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.MainPanel;
import org.broad.igv.ui.panel.Paintable;
import org.broad.igv.ui.svg.SVGGraphics;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.*;

/**
 * Utility methods for supporting saving of images as jpeg, png, and svg files.
 *
 * @author eflakes
 * @modified jrobinso
 */
public class SnapshotUtilities {

    /**
     * Class logger
     */
    private static Logger log = Logger.getLogger(SnapshotUtilities.class);


    /**
     * The maximum height in pixels for snapshots of a panel.
     */
    public static int DEFAULT_MAX_PANEL_HEIGHT = 1000;

    /**
     * We need to use a static for max panel height,  or alternatively much refactoring, however this class might
     * be accessed by multiple threads which set this to different values => use a thread local
     */
    private static ThreadLocal<Integer> maxPanelHeight = new ThreadLocal() {
        @Override
        protected Object initialValue() {
            return new Integer(DEFAULT_MAX_PANEL_HEIGHT);
        }
    };

    public static int getMaxPanelHeight() {
        return maxPanelHeight.get().intValue();
    }

    public static void setMaxPanelHeight(int h) {
        maxPanelHeight.set(h);
    }

    // Treat this class as a singleton, no instances allowed
    private SnapshotUtilities() {
    }


    public static void doComponentSnapshot(Component component, File file, SnapshotFileChooser.SnapshotFileType type) throws IOException{

        int width = component.getWidth();
        int height = component.getHeight();

        // Call appropriate converter
        switch (type) {
            case JPEG:
                exportScreenShotJPEG(component, file, width, height);
                break;
            //case EPS:
            //    exportScreenShotEPS(component, file, width, height);
            //    break;
            case PNG:
                exportScreenShotPNG(component, file, width, height);
                break;
            case SVG:
                log.debug("Exporting svg screenshot");
                exportScreenshotSVG(component, file);
                break;
        }
    }

    private static void exportScreenshotSVG2(Component target, File selecteddFile) {
        PrintWriter pw = null;
        try {
            pw = new PrintWriter(new FileWriter(selecteddFile));

            pw.println("<?xml version=\"1.0\" standalone=\"no\"?>\n" +
                    "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \n" +
                    "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n" +
                    "\n" +
                    "<svg width=\"100%\" height=\"100%\" version=\"1.1\"\n" +
                    "xmlns=\"http://www.w3.org/2000/svg\">");


            // TODO -- rectangle
            Rectangle rectangle = target.getBounds();
            SVGGraphics g2d = new SVGGraphics(pw, rectangle);

            target.paint(g2d);
            pw.print("</svg>");
            pw.close();
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }


    }


    private static void exportScreenshotSVG(Component target, File selectedFile) throws IOException {

        // Disable extending panel height beyond visible area
        //SnapshotUtilities.setMaxPanelHeight(-1);
        log.debug("Getting dom");
        DOMImplementation domImpl = GenericDOMImplementation.getDOMImplementation();


        // Create an instance of org.w3c.dom.Document.
        String svgNS = "http://www.w3.org/2000/svg";
        Document document = domImpl.createDocument(svgNS, "svg", null);

        // Create an instance of the SVG Generator.
        SVGGraphics2D svgGenerator = new SVGGraphics2D(document);
        //logger.info("Painting");
        target.paint(svgGenerator);

        Writer out = null;
        try {
            // Finally, stream out SVG to the standard output using
            // UTF-8 encoding.
            boolean useCSS = true; // we want to use CSS style attributes
            out = new BufferedWriter(new FileWriter(selectedFile));
            //logger.info("Writing output");
            svgGenerator.stream(out, useCSS);
        } finally {
            if (out != null) try {
                out.close();
            } catch (IOException e) {
                log.error("Error closing svg file", e);
            }
        }

    }

    private static void exportScreenShotJPEG(Component target, File selectedFile, int width, int height) throws IOException{

        BufferedImage image = getDeviceCompatibleImage(width, height);
        Graphics g = image.createGraphics();
        target.paintAll(g);

        String[] exts = new String[]{".jpg", ".jpeg"};
        boolean hasExt = false;
        if (selectedFile != null) {

            for(String ext: exts){
                if (selectedFile.getName().toLowerCase().endsWith(ext)) {
                    hasExt = true;
                    break;
                }
            }
            if(!hasExt){
                String correctedFilename = selectedFile.getAbsolutePath() + ".jpeg";
                selectedFile = new File(correctedFilename);
            }
            ImageIO.write(image, "jpeg", selectedFile);

        }
    }

    private static void exportScreenShotPNG(Component target, File selectedFile, int width, int height) throws IOException{

        BufferedImage image = getDeviceCompatibleImage(width, height);
        Graphics g = image.createGraphics();
        target.paintAll(g);

        if (selectedFile != null) {

            if (!selectedFile.getName().toLowerCase().endsWith(".png")) {
                String correctedFilename = selectedFile.getAbsolutePath() + ".png";
                selectedFile = new File(correctedFilename);
            }
            ImageIO.write(image, "png", selectedFile);


        }
    }


    /**
     * Creates a device compatible buffered svg.
     *
     * @param width  the svg width in pixels
     * @param height the svg height in pixels
     */
    public static BufferedImage getDeviceCompatibleImage(int width, int height) {

        GraphicsEnvironment graphicsEnvironment = GraphicsEnvironment.getLocalGraphicsEnvironment();
        GraphicsDevice screenDevice = graphicsEnvironment.getDefaultScreenDevice();
        GraphicsConfiguration graphicConfiguration = screenDevice.getDefaultConfiguration();
        BufferedImage image = graphicConfiguration.createCompatibleImage(width, height);

        return image;
    }


}
