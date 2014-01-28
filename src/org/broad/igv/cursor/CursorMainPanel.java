/*
 * Created by JFormDesigner on Tue Jan 14 11:39:00 EST 2014
 */

package org.broad.igv.cursor;

import java.awt.*;
import java.io.Serializable;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.plaf.*;

/**
 * @author Jim Robinson
 */
public class CursorMainPanel extends JPanel implements Serializable {

    private CursorModel model;

    public CursorMainPanel() {
        initComponents();
        //this.trackLabelPanel.setVmargin(10);
    }


    public void setModel(CursorModel model) {
        this.model = model;
        this.cursorIdeogramPanel1.setModel(model);
    }


    public void addTrack(CursorTrack track) {

        CursorTrackPanel panel = new CursorTrackPanel(track, model, this);
        this.trackPanel.add(panel);

        CursorTrackLabelPanel labelPanel = new CursorTrackLabelPanel(track, model, this);
        //JButton labelPanel = new JButton(track.getName());
        this.trackLabelPanel.add(labelPanel);

        this.cursorIdeogramPanel1.addTrack(track);
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        panel2 = new JPanel();
        cursorIdeogramPanel1 = new CursorIdeogramPanel();
        panel3 = new JLabel();
        scrollPane1 = new JScrollPane();
        panel1 = new JPanel();
        cursorAttributePanel1 = new JPanel();
        trackLabelPanel = new CursorTrackPanelContainer();
        cursorDataPanel1 = new JPanel();
        trackPanel = new CursorTrackPanelContainer();

        //======== this ========
        setBackground(Color.white);
        setLayout(new BorderLayout());

        //======== panel2 ========
        {
            panel2.setMinimumSize(new Dimension(24, 20));
            panel2.setLayout(new BorderLayout(4, 4));
            panel2.add(cursorIdeogramPanel1, BorderLayout.CENTER);

            //---- panel3 ----
            panel3.setPreferredSize(new Dimension(100, 40));
            panel2.add(panel3, BorderLayout.WEST);
        }
        add(panel2, BorderLayout.NORTH);

        //======== scrollPane1 ========
        {
            scrollPane1.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);

            //======== panel1 ========
            {
                panel1.setLayout(new BorderLayout(4, 4));

                //======== cursorAttributePanel1 ========
                {
                    cursorAttributePanel1.setBorder(LineBorder.createBlackLineBorder());
                    cursorAttributePanel1.setPreferredSize(new Dimension(100, 42));
                    cursorAttributePanel1.setLayout(new BorderLayout());
                    cursorAttributePanel1.add(trackLabelPanel, BorderLayout.CENTER);
                }
                panel1.add(cursorAttributePanel1, BorderLayout.WEST);

                //======== cursorDataPanel1 ========
                {
                    cursorDataPanel1.setLayout(new BorderLayout());

                    //---- trackPanel ----
                    trackPanel.setBorder(null);
                    trackPanel.setBackground(SystemColor.window);
                    cursorDataPanel1.add(trackPanel, BorderLayout.CENTER);
                }
                panel1.add(cursorDataPanel1, BorderLayout.CENTER);
            }
            scrollPane1.setViewportView(panel1);
        }
        add(scrollPane1, BorderLayout.CENTER);
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel panel2;
    private CursorIdeogramPanel cursorIdeogramPanel1;
    private JLabel panel3;
    private JScrollPane scrollPane1;
    private JPanel panel1;
    private JPanel cursorAttributePanel1;
    private CursorTrackPanelContainer trackLabelPanel;
    private JPanel cursorDataPanel1;
    private CursorTrackPanelContainer trackPanel;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
