package org.igv.ui;

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
