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

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.util.MessageUtils;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

import java.awt.*;
import java.awt.geom.AffineTransform;
import java.io.*;

/**
 * @author jrobinso
 * @date Nov 20, 2010
 */
public class SVGTest {

    public static void main(String[] args) throws IOException {

        test2();


    }

    private static void test1() throws IOException {
        PrintWriter pw = new PrintWriter(new FileWriter("test.svg"));

        pw.println("<?xml version=\"1.0\" standalone=\"no\"?>\n" +
                "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \n" +
                "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n" +
                "\n" +
                "<svg width=\"100%\" height=\"100%\" version=\"1.1\"\n" +
                "xmlns=\"http://www.w3.org/2000/svg\">");


        SVGGraphics g2d = new SVGGraphics(pw, new Rectangle(0,0,500,500));

        draw(g2d);

        pw.print("</svg>");

        pw.close();
    }

    private static void draw(Graphics2D g2d) {
        g2d.clipRect(5, 5, 350, 350);

        g2d.setFont(FontManager.getFont(Font.BOLD, 24));
        g2d.drawString("Hello World", 0, 20);


        g2d.setColor(Color.red);
        g2d.fillRect(10, 10, 50, 100);

        g2d.drawRect(100, 10, 50, 100);

        g2d.setColor(Color.blue);
        g2d.drawRect(100, 60, 100, 50);

        g2d.setColor(Color.lightGray);
        int[] x = {220, 300, 170, 123};
        int[] y = {100, 210, 250, 234};
        g2d.fillPolygon(x, y, x.length);


        g2d.setColor(Color.green);
        g2d.drawLine(200, 200, 400, 400);

        AffineTransform tr = new AffineTransform();
        tr.rotate(Math.PI / 2);
        tr.translate(0, -20);
        g2d.setTransform(tr);


        g2d.setColor(Color.black);
        g2d.drawString("Hello World", 0, 0);
    }

    private static void test2() {

        File selecteddFile = new File("test2.svg");
         try {

             DOMImplementation domImpl = GenericDOMImplementation.getDOMImplementation();


             // Create an instance of org.w3c.dom.Document.
             String svgNS = "http://www.w3.org/2000/svg";
             Document document = domImpl.createDocument(svgNS, "svg", null);


             // Create an instance of the SVG Generator.
             org.apache.batik.svggen.SVGGraphics2D svgGenerator = new SVGGraphics2D(document);

             draw(svgGenerator);

             Writer out = new BufferedWriter(new FileWriter(selecteddFile));
             //logger.info("Writing output");
             svgGenerator.stream(out, false);
             //logger.info("Done");
         } catch (Exception e) {
             MessageUtils.showMessage("Error encountered creating SVG file: " + e.toString());
         }
     }

}
