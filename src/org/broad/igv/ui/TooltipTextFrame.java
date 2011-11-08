package org.broad.igv.ui;

import org.broad.igv.util.BrowserLauncher;

import java.awt.*;
import java.awt.event.*;
import java.io.IOException;
import javax.swing.*;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;

/**
 * @author Jim Robinson
 * @date 11/7/11
 */
public class TooltipTextFrame extends JFrame {


    public static void main(String[] args) {

        StringBuffer text = new StringBuffer("<html><b><i>Hello Word</b></i><br> " +
                "<a href=\"http://www.google.com\">Google</a><br>");
        for (int i = 0; i < 100; i++) {
            for (int j = 0; j < 40; j++) {
                text.append("blah ");
            }
            text.append(".....");
            text.append("<br>");
        }
        text.append("</html>");

        TooltipTextFrame frame = new TooltipTextFrame(text.toString());
        frame.setVisible(true);

    }

    public TooltipTextFrame(String text) throws HeadlessException {
        init(text);
    }

    void init(String text) {

        setUndecorated(true);
        setAlwaysOnTop(true);

        JButton button = new JButton("Close Me");
        button.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                System.exit(0);
            }
        });


        setSize(300, 300);
        setLocation(200, 200);
        setLayout(new BorderLayout());

        getContentPane().add(button, BorderLayout.NORTH);


        JEditorPane pane = new JEditorPane("text/html", text);
        pane.setEditable(false);
        pane.setBackground(new Color(255, 254, 217));
        JScrollPane scrollPane = new JScrollPane(pane);

        MouseAdapter mouseAdapter = new MouseAdapter() {

            private Point point = new Point();

            @Override
            public void mouseDragged(MouseEvent e) {
                Point p = getLocation();
                final int dx = e.getX() - point.x;
                final int dy = e.getY() - point.y;
                setLocation(p.x + dx, p.y + dy);
            }

            @Override
            public void mousePressed(MouseEvent e) {
                point.x = e.getX();
                point.y = e.getY();
            }
        };

        pane.addMouseListener(mouseAdapter);

        pane.addMouseMotionListener(mouseAdapter);

        pane.addHyperlinkListener(new HyperlinkListener() {
            public void hyperlinkUpdate(HyperlinkEvent e) {
                try {
                    if (e.getEventType() == HyperlinkEvent.EventType.ACTIVATED)
                        BrowserLauncher.openURL(e.getURL().toExternalForm());
                } catch (IOException e1) {
                    e1.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }
        });

        getContentPane().add(scrollPane, BorderLayout.CENTER);
    }


}
