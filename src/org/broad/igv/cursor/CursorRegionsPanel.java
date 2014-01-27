package org.broad.igv.cursor;

import org.broad.igv.feature.*;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.util.Pair;
import org.broad.tribble.Feature;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * @author jrobinso
 *         Date: 1/14/14
 *         Time: 12:41 PM
 */
public class CursorRegionsPanel extends JComponent implements Serializable {

    Cursor model;
    Color regionGray = new Color(193, 193, 193);
    Color lightBlue = new Color(0, 0, 150);
    private CursorMainPanel mainPanel;

    public CursorRegionsPanel() {
        setBorder(BorderFactory.createLineBorder(Color.black));

        MouseAdapter ma = new MouseHandler();
        addMouseListener(ma);
        addMouseMotionListener(ma);
    }


    public void setModel(Cursor model) {
        this.model = model;
    }

    @Override
    protected void paintComponent(Graphics graphics) {
        super.paintComponent(graphics);

        if (model == null) return;

        List<CursorFrame> frameList = model.getFrames();
        int frameMargin = model.getFrameMargin();

        if(frameList == null) return;

        int framePixelWidth = model.getFramePixelWidth();
        int frameBPWidth = model.getFrameBPWidth();
        double origin = model.getOrigin();   // In frame units
        double end = origin + ((double) getWidth()) / framePixelWidth;

        int h = model.getTrackPixelHeight();


        int startFrameNumber = (int) origin;
        if (startFrameNumber >= frameList.size()) return;

        for (int frameNumber = startFrameNumber; frameNumber < frameList.size(); frameNumber++) {

            if (frameNumber > end) break;

            CursorFrame frame = frameList.get(frameNumber);
            double pxStart = (frameNumber - origin) * framePixelWidth;
            double pxEnd = pxStart + framePixelWidth;

            // Region block
            graphics.setColor(regionGray);
            graphics.fillRect((int) (pxStart + frameMargin/2),  0, framePixelWidth - frameMargin, h);
            graphics.setColor(lightBlue);
            int lineX = (int) ((pxStart + pxEnd) / 2);
            graphics.drawLine(lineX, 1, lineX, h-2);


        }
    }

    public void setMainPanel(CursorMainPanel mainPanel) {
        this.mainPanel = mainPanel;
    }

    public CursorMainPanel getMainPanel() {
        return mainPanel;
    }

    class MouseHandler extends MouseAdapter {

        int mouseX;

        @Override
        public void mouseClicked(MouseEvent mouseEvent) {

        }

        @Override
        public void mousePressed(MouseEvent mouseEvent) {
            mouseX = mouseEvent.getX();
        }

        @Override
        public void mouseDragged(MouseEvent mouseEvent) {
            int delta = mouseX - mouseEvent.getX();
            model.shiftOriginPixels(delta);
            mouseX = mouseEvent.getX();

            // TODO -- use event here and remove back pointer

            mainPanel.repaint();

        }
    }

}
