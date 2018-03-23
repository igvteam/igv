/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
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

        test1();


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
