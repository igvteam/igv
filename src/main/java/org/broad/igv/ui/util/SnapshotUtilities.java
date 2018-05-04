/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2018 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/**
 * SnapshotUtilities.java
 * <p>
 * Created on November 29, 2007, 2:14 PM
 * <p>
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package org.broad.igv.ui.util;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.apache.log4j.Logger;
import org.broad.igv.ui.panel.MainPanel;
import org.broad.igv.ui.panel.Paintable;
import org.broad.igv.util.RuntimeUtils;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.*;
import java.lang.reflect.Constructor;

/**
 * Utility methods for supporting saving of images as jpeg, png, and svg files.
 *
 * @author eflakes
 * @modified jrobinso
 */
public class SnapshotUtilities {

    private static Logger log = Logger.getLogger(SnapshotUtilities.class);


    private static String EPSClassName = "net.sf.epsgraphics.EpsGraphics";
    private static String EPSColorModeClassName = "net.sf.epsgraphics.ColorMode";


    /**
     * The maximum height in pixels for snapshots of a panel.
     */
    public static int DEFAULT_MAX_PANEL_HEIGHT = 1000;

    /**
     * We need to use a static for max panel height,  or alternatively much refactoring
     */
    private static int maxPanelHeight = DEFAULT_MAX_PANEL_HEIGHT;

    public static int getMaxPanelHeight() {
        return maxPanelHeight;
    }

    public static void setMaxPanelHeight(int h) {
        maxPanelHeight = h;
    }

    // Treat this class as a singleton, no instances allowed
    private SnapshotUtilities() {
    }


    public static String doComponentSnapshot(Component component, File file, SnapshotFileChooser.SnapshotFileType type, boolean paintOffscreen) throws IOException {

        //TODO Should really make this work for more components
        if (paintOffscreen && !(component instanceof Paintable)) {
            log.error("Component cannot be painted offscreen. Performing onscreen paint");
            paintOffscreen = false;
        }

        if (paintOffscreen) {

            Rectangle rect = component.getBounds();

            if (component instanceof MainPanel) {
                rect.height = ((MainPanel) component).getOffscreenImageHeight();
            } else {
                rect.height = Math.min(component.getHeight(), getMaxPanelHeight());
            }

            // translate to (0, 0) if necessary
            int dx = rect.x;
            int dy = rect.y;
            rect.x = 0;
            rect.y = 0;
            rect.width -= dx;
            rect.height -= dy;

            component.setBounds(rect);
        }

        int width = component.getWidth();
        int height = component.getHeight();

        // Call appropriate converter
        String format = null;
        String[] exts = null;
        switch (type) {
            case SVG:
                //log.debug("Exporting svg screenshot");
                exportScreenshotSVG(component, file, width, height, paintOffscreen);
                //exportScreenshotVector2D(component, file, paintOffscreen);
                break;
            case JPEG:
                format = "jpeg";
                exts = new String[]{".jpg", ".jpeg"};
                break;
            case PNG:
                format = "png";
                exts = new String[]{"." + format};
                break;
            case EPS:
                exportScreenshotEpsGraphics(component, file, width, height, paintOffscreen);
                //exportScreenshotEpsGraphicsNoRef(component, file, paintOffscreen);
                break;
        }
        if (format != null && exts != null) {
            exportScreenShotBufferedImage(component, file, width, height, exts, format, paintOffscreen);
        }
        return "OK";
    }

//    private static void exportScreenshotVector2D(Component target, File selectedFile, boolean paintOffscreen) throws IOException{
//
//        de.erichseifert.vectorgraphics2d.VectorGraphics2D g = null;
//        String filePath = selectedFile.getAbsolutePath();
//
//        if(filePath.endsWith(".svg")){
//            g = new de.erichseifert.vectorgraphics2d.SVGGraphics2D(0.0, 0.0, target.getWidth(), target.getHeight());
//        }else if(filePath.endsWith(".eps")){
//            g = new de.erichseifert.vectorgraphics2d.EPSGraphics2D(0.0, 0.0, target.getWidth(), target.getHeight());
//        }
//        target.paintAll(g);
//
//        // Write the output to a file
//        FileOutputStream file = new FileOutputStream(selectedFile);
//        try {
//            file.write(g.getBytes());
//        } finally {
//            file.close();
//        }
//
//    }

    private static void exportScreenshotEpsGraphics(Component target, File selectedFile, int width, int height, boolean paintOffscreen) throws IOException {

        if (!SnapshotUtilities.canExportScreenshotEps()) {
            String msg = "ERROR: EPS output requires EPSGraphics library. See https://www.broadinstitute.org/software/igv/third_party_tools#epsgraphics";
            log.error(msg);
            return;
        }

        Graphics2D g = null;
        FileOutputStream fos = null;
        try {
            Class colorModeClass = RuntimeUtils.loadClassForName(EPSColorModeClassName, null);
            Class graphicsClass = RuntimeUtils.loadClassForName(EPSClassName, null);

            Constructor constructor = graphicsClass.getConstructor(String.class, OutputStream.class,
                    int.class, int.class, int.class, int.class, colorModeClass);

            Object colorModeValue = Enum.valueOf(colorModeClass, "COLOR_RGB");

            // EpsGraphics stores directly in a file
            fos = new FileOutputStream(selectedFile);
            g = (Graphics2D) constructor.newInstance("eps", fos, 0, 0, target.getWidth(), target.getHeight(), colorModeValue);

            paintImage(target, g, width, height, paintOffscreen);

            graphicsClass.getMethod("close").invoke(g);

        } catch (Exception e) {
            log.error(e.getMessage(), e);
        } finally {
            if (fos != null) {
                fos.flush();
                fos.close();
            }
        }

    }

    public static boolean canExportScreenshotEps() {
        Constructor constr = null;
        try {
            Class colorModeClass = RuntimeUtils.loadClassForName(EPSColorModeClassName, null);
            Class graphicsClass = RuntimeUtils.loadClassForName(EPSClassName, null);
            constr = graphicsClass.getConstructor(String.class, OutputStream.class,
                    int.class, int.class, int.class, int.class, colorModeClass);
        } catch (Exception e) {
            //pass
        }
        return constr != null;
    }

//    private static void exportScreenshotEpsGraphicsNoRef(Component target, File selectedFile, boolean paintOffscreen) throws IOException{
//
//        FileOutputStream fos = null;
//        try {
//            // EpsGraphics stores directly in a file
//            fos = new FileOutputStream(selectedFile);
//            Graphics2D g = new EpsGraphics("eps", fos, 0, 0, target.getWidth(), target.getHeight(), ColorMode.COLOR_RGB);
//            paintImage(target, g, paintOffscreen);
//
//        } finally {
//            if(fos != null){
//                fos.flush();
//                fos.close();
//            }
//        }
//
//    }

    private static void exportScreenshotSVG(Component target, File selectedFile, int width, int height, boolean paintOffscreen) throws IOException {

        String format = "svg";
        selectedFile = fixFileExt(selectedFile, new String[]{format}, format);

        // Create an instance of org.w3c.dom.Document.                                                                                      
        DOMImplementation domImpl = GenericDOMImplementation.getDOMImplementation();
        String svgNS = "http://www.w3.org/2000/svg";
        Document document = domImpl.createDocument(svgNS, format, null);

        // Write image data into document                                                                                                   
        SVGGraphics2D svgGenerator = new SVGGraphics2D(document);

        paintImage(target, svgGenerator, width, height, paintOffscreen);

        Writer out = null;
        try {
            // Finally, stream out SVG to the standard output using                                                                         
            // UTF-8 encoding.                                                                                                              
            boolean useCSS = true; // we want to use CSS style attributes                                                                   
            out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(selectedFile), "UTF-8"));
            svgGenerator.stream(out, useCSS);
        } finally {
            if (out != null) try {
                out.close();
            } catch (IOException e) {
                log.error("Error closing svg file", e);
            }
        }
    }

    private static void paintImage(Component target, Graphics2D g, int width, int height, boolean paintOffscreen) {
        log.debug("Painting to target " + target + " , offscreen " + paintOffscreen);
        if (paintOffscreen) {
            Rectangle rect = new Rectangle(0, 0, width, height);
            ((Paintable) target).paintOffscreen(g, rect);
        } else {
            target.paintAll(g);
        }
    }

    /**
     * Export the specified {@code target} component as a {@code BufferedImage} to the given file.
     *
     * @param target
     * @param selectedFile
     * @param width
     * @param height
     * @param allowedExts
     * @param format         Format, also appended as an extension if the file doesn't end with anything in {@code allowedExts}
     * @param paintOffscreen
     * @throws IOException
     */
    private static void exportScreenShotBufferedImage(Component target, File selectedFile, int width, int height,
                                                      String[] allowedExts, String format, boolean paintOffscreen) throws IOException {
        BufferedImage image = getDeviceCompatibleImage(width, height);
        Graphics2D g = image.createGraphics();

        paintImage(target, g, width, height, paintOffscreen);

        selectedFile = fixFileExt(selectedFile, allowedExts, format);
        if (selectedFile != null) {
            log.debug("Writing image to " + selectedFile.getAbsolutePath());
            ImageIO.write(image, format, selectedFile);
        }
    }

    public static BufferedImage createBufferedImage(Component target, Rectangle clipRect,  int maxHeight) throws IOException {

        int width = clipRect.width;
        int height = Math.min(clipRect.height, 1000);

        BufferedImage image = getDeviceCompatibleImage(width, height);
        Graphics2D g = image.createGraphics();

       // clipRect.x += 5;
       // clipRect.width -=5;

       // g.setClip(clipRect);
        g.translate(0, -clipRect.y);
        height += clipRect.y;

        Rectangle rect = new Rectangle(0, 0, width, height);
        ((Paintable) target).paintOffscreen(g, rect);

        g.setBackground(Color.WHITE);
        g.clearRect(0, 0, 5, height);

        return image;

    }

    /**
     * Add a file extension to the file if it doesn't already
     * have an acceptable one
     *
     * @param selectedFile
     * @param allowedExts  Strings which qualify as extensions
     * @param defExtension Default extension. A period be inserted in between the file path iff {@code defExtension}
     *                     does not already have it
     * @return Either the input File, if it had an extension contained in {@code allowedExts},
     * or a new with with {@code defExtension} appended
     */
    private static File fixFileExt(File selectedFile, String[] allowedExts, String defExtension) {
        boolean hasExt = false;
        if (selectedFile != null) {
            for (String ext : allowedExts) {
                if (selectedFile.getName().toLowerCase().endsWith(ext)) {
                    hasExt = true;
                    break;
                }
            }
            if (!hasExt) {
                String addExt = defExtension.startsWith(".") ? defExtension : "." + defExtension;
                String correctedFilename = selectedFile.getAbsolutePath() + addExt;
                selectedFile = new File(correctedFilename);
            }
        }
        return selectedFile;
    }


    /**
     * Creates a device compatible BufferedImage
     *
     * @param width  the width in pixels
     * @param height the height in pixels
     */
    public static BufferedImage getDeviceCompatibleImage(int width, int height) {

        GraphicsEnvironment graphicsEnvironment = GraphicsEnvironment.getLocalGraphicsEnvironment();
        GraphicsDevice screenDevice = graphicsEnvironment.getDefaultScreenDevice();
        GraphicsConfiguration graphicConfiguration = screenDevice.getDefaultConfiguration();
        BufferedImage image = graphicConfiguration.createCompatibleImage(width, height);

        return image;
    }


}
