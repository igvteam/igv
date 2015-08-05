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

package org.broad.igv.synteny;

import org.broad.igv.util.ChromosomeColors;

import java.awt.Color;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;

public class Region extends AbstractMapping
{
  List<Anchor> anchors = new ArrayList();

  public void addAnchor(Anchor a)
  {
    String fromDirection = "+";
    String fromChr = getFromChr();
    String toDirection = a.getDirection() ? "+" : "-";
    String toChr = a.getToChr();

    if (this.anchors.size() > 0) {
      Anchor lastAnchor = anchors.get(this.anchors.size() - 1);
      int lastFromEnd = lastAnchor.getFromEnd();
      int lastToEnd = lastAnchor.getToEnd();

      Anchor fillerAnchor = new Anchor();
      int fromStart = lastFromEnd;
      int fromEnd = a.getFromStart();
      int toStart = lastToEnd;
      int toEnd = a.getToStart();
      if ((toEnd > toStart) && (fromEnd > fromStart)) {
        fillerAnchor.setParameters("psueudo", fromChr, fromStart, fromEnd, fromDirection,
                toChr, lastToEnd, a.getToStart(), toDirection);
        fillerAnchor.psuedo = true;
        this.anchors.add(fillerAnchor);
      }
    }

    this.anchors.add(a);
  }

  public List<Anchor> getAnchors() {
    return this.anchors;
  }

  public double mapPosition(int position)
  {
    for (Anchor a : this.anchors) {
      if (a.containsFromPosition(position)) {
        return a.mapPosition(position);
      }
    }

    if (containsFromPosition(position)) {
      double delta = this.scaleFactor * (position - this.fromStart);

      if (this.direction == true) {
        return this.toStart + delta;
      }
      return this.toEnd - delta;
    }

    return -1.0D;
  }

  public String toBed()
  {
    StringWriter sw = new StringWriter();
    PrintWriter pw = new PrintWriter(sw);
    String name = getToChr() + ":" + getToStart() + "-" + getToEnd();

    int blockCount = 0;
    String blockSizes = "";
    String blockStarts = "";
    for (Anchor a : this.anchors) {
      if (!a.psuedo) {
        blockCount++;
        blockSizes = blockSizes + (a.getFromEnd() - a.getFromStart()) + ",";
        blockStarts = blockStarts + (a.getFromStart() - getFromStart()) + ",";
      }
    }

    pw.print(getFromChr() + "\t" + getFromStart() + "\t" + getFromEnd() + "\t" + name + "\t1000\t" +
            (getDirection() ? "+" : "-") + "\t" + getFromStart() + "\t" + getFromEnd() + "\t" +
            convertColorToRGBString(ChromosomeColors.getColor(getToChr())) + "\t" + blockCount + "\t" +
            blockSizes + "\t" + blockStarts);

    return sw.toString();
  }

  public static String convertColorToRGBString(Color color)
  {
    StringBuffer buffer = new StringBuffer();
    buffer.append(color.getRed());
    buffer.append(",");
    buffer.append(color.getGreen());
    buffer.append(",");
    buffer.append(color.getBlue());
    return buffer.toString();
  }
}