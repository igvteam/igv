package org.broad.igv.cursor;


import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.Serializable;

/**
 * @author jrobinso
 *         Date: 1/26/14
 *         Time: 10:55 PM
 */
public class CursorTrackLabelPanel extends JComponent implements Serializable {

    final CursorModel model;
    final CursorTrack track;
    int sortDirection = 1;

    JButton sortButton;
    CursorIdeogramPanel ideogramPanel;

    public CursorTrackLabelPanel(final CursorTrack track, final CursorModel model, final CursorMainPanel mainPanel) {

        this.track = track;
        this.model = model;
        this.sortButton = new JButton("Sort");
        this.sortButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                model.sortFrames(track, sortDirection);
                sortDirection *= -1;
                mainPanel.repaint();
            }
        });
        add(sortButton);

        ideogramPanel = new CursorIdeogramPanel();
        ideogramPanel.addTrack(track);
        ideogramPanel.setModel(model);
        ideogramPanel.setDrawViewRect(false);
        add(ideogramPanel);
    }


    @Override
    public void doLayout() {
        synchronized (this.getTreeLock()) {
            Rectangle bound = getBounds();
            sortButton.setBounds(0, 0, 50, bound.height);
            ideogramPanel.setBounds(50, 0, bound.width - 50, bound.height);
        }
    }
}
