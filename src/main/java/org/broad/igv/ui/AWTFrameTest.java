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
