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

    ButtonGroup trackSelectionButtonGroup;


    public CursorMainPanel() {
        initComponents();
        trackSelectionButtonGroup = new ButtonGroup();
        //this.trackLabelPanel.setVmargin(10);
    }


    public void addTrackSelectionButton(JRadioButton button) {
        if(trackSelectionButtonGroup.getButtonCount() == 0) {
            button.setSelected(true);
        }
        trackSelectionButtonGroup.add(button);
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

    @Override
    public void revalidate() {
        super.revalidate();
    }

    public void tracksAdded() {
        if(this.scrollPane1 != null) this.scrollPane1.revalidate();
    }

    public int getDataPanelWidth() {
        return dataPanel.getWidth();
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        panel2 = new JPanel();
        cursorIdeogramPanel1 = new CursorIdeogramPanel();
        panel3 = new JLabel();
        scrollPane1 = new JScrollPane();
        panel1 = new JPanel();
        attributePanel = new JPanel();
        trackLabelPanel = new CursorTrackPanelContainer();
        dataPanel = new JPanel();
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

                //======== attributePanel ========
                {
                    attributePanel.setBorder(LineBorder.createGrayLineBorder());
                    attributePanel.setPreferredSize(new Dimension(100, 42));
                    attributePanel.setLayout(new BorderLayout());
                    attributePanel.add(trackLabelPanel, BorderLayout.CENTER);
                }
                panel1.add(attributePanel, BorderLayout.WEST);

                //======== dataPanel ========
                {
                    dataPanel.setBorder(LineBorder.createGrayLineBorder());
                    dataPanel.setLayout(new BorderLayout());

                    //---- trackPanel ----
                    trackPanel.setBorder(null);
                    trackPanel.setBackground(SystemColor.window);
                    dataPanel.add(trackPanel, BorderLayout.NORTH);
                }
                panel1.add(dataPanel, BorderLayout.CENTER);
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
    private JPanel attributePanel;
    private CursorTrackPanelContainer trackLabelPanel;
    private JPanel dataPanel;
    private CursorTrackPanelContainer trackPanel;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
