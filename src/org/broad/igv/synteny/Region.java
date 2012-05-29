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