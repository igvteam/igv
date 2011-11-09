package org.broad.igv.ui;

import org.broad.igv.util.BrowserLauncher;
import org.broad.igv.util.HttpUtils;

import java.awt.*;
import java.awt.event.*;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.net.MalformedURLException;
import java.net.URL;
import javax.swing.*;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;

/**
 * @author Jim Robinson
 * @date 11/7/11
 */
public class TooltipTextFrame extends JFrame {


    public static void main(String[] args) throws IOException {

        test3();
        //omaTest();

    }

    private static void test3() {

        JFrame frame = new JFrame();
        frame.setSize(300, 300);
        JLabel label = new JLabel("Hello");
        label.setToolTipText("<html><b><i>Hello Word</b></i><br> " +
                "<a href=\"http://www.google.com\">Google</a><br>\n");
        frame.add(label);
        frame.setVisible(true);

    }

    private static void test1() {
        StringBuffer text = new StringBuffer("<html><b><i>Hello Word</b></i><br> " +
                "<a href=\"http://www.google.com\">Google</a><br>\n");
        for (int i = 0; i < 100; i++) {
            for (int j = 0; j < 40; j++) {
                text.append("blah ");
            }
            text.append("<div style=\"color:red\"> .....</div>");
            text.append("<br>\n");
        }
        text.append("</html>");

        TooltipTextFrame frame = new TooltipTextFrame(text.toString());
        frame.setVisible(true);

    }

    public static void omaTest() throws IOException {
        String url = "http://mutationassessor.org/v1/?cm=var&var=hg18,7,55178574,G,A&frm=txt";
        String result = HttpUtils.getInstance().getContentsAsString(new URL(url));

        BufferedReader br = new BufferedReader(new StringReader(result));
        String[] headers = br.readLine().split("\t");
        String[] values = br.readLine().split("\t");

        StringBuffer buf = new StringBuffer();
        buf.append("<html>");

        buf.append("<table>");
        int n = Math.min(headers.length, values.length);
        for (int i = 0; i < n; i++) {
            buf.append("<tr>");
            buf.append("<td>");
            final String header = headers[i];
            buf.append(header);
            buf.append("</td>");

            if (header.equals("MSA") || header.equals("PDB")) {
                buf.append("<a href=\"http://" + values[i] + "\">");
                buf.append(header);
                buf.append("</a>");
            } else if (header.equals("Uniprot")) {
                buf.append("<a href=\"http://www.uniprot.org/uniprot/" + values[i] + "\">");
                buf.append(values[i]);
                buf.append("</a>");
            } else if (header.equals("Refseq")) {
                buf.append("<a href=\"http://www.ncbi.nlm.nih.gov/sites/entrez?db=protein&cmd=search&term=" + values[i] + "\">");
                buf.append(values[i]);
                buf.append("</a>");
            } else {
                buf.append(values[i]);
            }
            buf.append("</td></tr>");
        }
        buf.append("</table>");


        TooltipTextFrame frame = new TooltipTextFrame(buf.toString());
        frame.setVisible(true);

    }

    public TooltipTextFrame(String text) throws HeadlessException {
        init(text);
    }

    void init(String text) {

        //setUndecorated(true);
        setAlwaysOnTop(true);

        setSize(300, 300);
        setLayout(new BorderLayout());

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
