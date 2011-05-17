/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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

/*
 * Created by JFormDesigner on Tue May 17 18:20:38 EDT 2011
 */

package org.broad.igv.peaks;

import org.broad.igv.ui.IGV;

import java.awt.*;
import javax.swing.*;
import javax.swing.event.*;

/**
 * @author Stan Diamond
 */
public class PeakCommandBar extends JPanel {
    public PeakCommandBar() {
        initComponents();
        scoreSlider.setValue((int) PeakTrack.getScoreThreshold());
        foldChangeSlider.setValue((int) PeakTrack.getFoldChangeThreshold());
    }

    private void foldChangeSliderStateChanged(ChangeEvent e) {
        PeakTrack.setScoreThreshold(foldChangeSlider.getValue());
        IGV.getInstance().repaintDataPanels();

    }

    private void scoreSliderStateChanged(ChangeEvent e) {
        PeakTrack.setFoldChangeThreshold(scoreSlider.getValue());
        IGV.getInstance().repaintDataPanels();
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        label1 = new JLabel();
        foldChangeSlider = new JSlider();
        label2 = new JLabel();
        scoreSlider = new JSlider();

        //======== this ========
        setLayout(new FlowLayout(FlowLayout.LEFT));

        //---- label1 ----
        label1.setText("Fold change");
        add(label1);

        //---- foldChangeSlider ----
        foldChangeSlider.setPaintTicks(true);
        foldChangeSlider.setPaintLabels(true);
        foldChangeSlider.setToolTipText("Adjust score threshold");
        foldChangeSlider.setMajorTickSpacing(2);
        foldChangeSlider.setMinorTickSpacing(1);
        foldChangeSlider.setMaximum(10);
        foldChangeSlider.setValue(0);
        foldChangeSlider.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                foldChangeSliderStateChanged(e);
            }
        });
        add(foldChangeSlider);

        //---- label2 ----
        label2.setText("Score");
        add(label2);

        //---- scoreSlider ----
        scoreSlider.setPaintTicks(true);
        scoreSlider.setPaintLabels(true);
        scoreSlider.setToolTipText("Adjust score threshold");
        scoreSlider.setMajorTickSpacing(20);
        scoreSlider.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                scoreSliderStateChanged(e);
            }
        });
        add(scoreSlider);
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JLabel label1;
    private JSlider foldChangeSlider;
    private JLabel label2;
    private JSlider scoreSlider;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
