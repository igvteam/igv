package org.broad.igv.ui.panel;

import org.broad.igv.lists.GeneList;
import org.broad.igv.ui.IGV;

import javax.swing.*;
import javax.swing.border.Border;
import java.awt.*;

/**
 * @author Jim Robinson
 * @date 5/25/12
 */
public class PanelLayoutManager {

//    HeaderPanelContainer headerPanelContainer;
//    DataPanelContainer dataPanelContainer;
//
//    public void layoutContainer(Container container) {
//        synchronized (container.getTreeLock()) {
//            Component[] children = container.getComponents();
//            //java.util.List<ReferenceFrame> frames = FrameManager.getFrames();
//            int h = container.getHeight();
//
//            try {
//                //Sometimes the number of children is not the same as the number of frames.
//                //Not entirely sure why - Jacob S
//                for (int i = 0; i < Math.min(children.length, frames.size()); i++) {
//
//                    Component c = children[i];
//                    ReferenceFrame frame = frames.get(i);
//                    c.setBounds(frame.pixelX, 0, frame.getWidthInPixels(), h);
//
//                    if (c instanceof JComponent) {
//                        if (frame.getWidthInPixels() > 5) {
//                            ((JComponent) c).setBorder(panelBorder);
//                        } else {
//                            ((JComponent) c).setBorder(null);
//                        }
//                    }
//
//                    log.debug("Layout: " + frame.getName() + "  x=" + frame.pixelX + "  w=" + frame.getWidthInPixels());
//
//                }
//            } catch (Exception e) {
//                log.error("Error laying out data panel", e);
//            }
//        }
//    }
//
//    public void createHeaderPanels() {
//
//        removeAll();
//        frames = FrameManager.getFrames();
//
//        contentPanel = new JPanel();
//        contentPanel.setLayout(new DataPanelLayout());
//        for (ReferenceFrame f : frames) {
//            HeaderPanel dp = new HeaderPanel(f);
//            dp.setBackground(getBackground());
//            contentPanel.add(dp);
//        }
//        add(contentPanel, BorderLayout.CENTER);
//
//        if (FrameManager.isGeneListMode()) {
//            GeneList gl = IGV.getInstance().getSession().getCurrentGeneList();
//            String name = gl.getDisplayName();
//            JLabel label = new JLabel(name, JLabel.CENTER);
//            Border border = BorderFactory.createLineBorder(Color.lightGray);
//            label.setBorder(border);
//            add(label, BorderLayout.NORTH);
//        }
//
//        invalidate();
//    }
//
//    public void createDataPanels() {
//        removeAll();
//        for (ReferenceFrame f : FrameManager.getFrames()) {
//            DataPanel dp = new DataPanel(f, this);
//            add(dp);
//        }
//        invalidate();
//    }

}
