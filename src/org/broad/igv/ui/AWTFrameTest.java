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

package org.broad.igv.ui;

import javax.swing.*;
import java.awt.*;

/**
 * @author jrobinso
 * @date Mar 31, 2011
 */
public class AWTFrameTest {

    public static void main(String[] args) {



        JFrame jf = new JFrame();
        Component c = jf.getContentPane();  // <= a JPanel with a JRootPane$1 layout manager.
        Component c2 = jf.getRootPane().getContentPane();
        boolean isEqual = (c == c2);

        JRootPane jrp = new JRootPane();
        
        LayoutManager lm = jrp.getLayout();

        LayoutManager lm1 = jf.getContentPane().getLayout(); // <= a hcak to get a JRootPane$1 instance.  Maybe not neccessary?

        final Frame frame = new Frame();
        frame.setSize(400, 300);

        final JRootPane rootPane = new JRootPane();
        rootPane.setSize(frame.getSize());
        //rootPanel.setBackground(Color.blue);

        rootPane.setLayout(new BorderLayout());
        rootPane.add(new TestPanel());

        frame.add(rootPane);
       

        

        frame.setVisible(true);


    }

    static class TestPanel extends JPanel {

        @Override
        protected void paintComponent(Graphics graphics) {
            graphics.setColor(Color.cyan);
            graphics.fillRect(0, 0, getWidth(), getHeight()) ;
        }
    }


}
