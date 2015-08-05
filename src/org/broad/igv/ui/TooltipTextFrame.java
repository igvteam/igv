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

import org.broad.igv.util.BrowserLauncher;
import org.broad.igv.util.HttpUtils;

import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.awt.event.*;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import javax.swing.*;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;
import javax.swing.text.MutableAttributeSet;
import javax.swing.text.html.HTML;
import javax.swing.text.html.HTMLEditorKit;
import javax.swing.text.html.parser.ParserDelegator;

/**
 * @author Jim Robinson
 * @date 11/7/11
 */
public class TooltipTextFrame extends JFrame {


    private static final DataFlavor[] supportedFlavors;

    static {
        try {
            supportedFlavors = new DataFlavor[]{
                    new DataFlavor("text/html;class=java.lang.String"),
                    new DataFlavor("text/plain;class=java.lang.String")
            };
        } catch (ClassNotFoundException e) {
            throw new ExceptionInInitializerError(e);
        }
    }


    public TooltipTextFrame(String title, String text) throws HeadlessException {

        setTitle(title);

        //setUndecorated(true);
        setAlwaysOnTop(true);


        setLayout(new BorderLayout());

        JEditorPane pane = new JEditorPane("text/html", text);
        pane.setEditable(false);
        pane.setTransferHandler(new MyTransferHandler());

        Dimension d = pane.getPreferredSize();
        int w = (int) (1.1 * d.width);
        int h = (int) (1.1 * d.height);

        h = h > 600 ? 600 : (h < 100 ? 100 : h);
        w = w > 800 ? 800 : (w < 100 ? 100 : w);
        setSize(w, h);


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


    /**
     * Custom transfer handler that preserves line breaks in html when copying to clipboard
     */
    class MyTransferHandler extends TransferHandler {

        protected Transferable createTransferable(JComponent c) {
            final JEditorPane pane = (JEditorPane) c;
            final String htmlText = pane.getText();
            final String plainText = extractText(new StringReader(htmlText));
            return new MyTransferable(plainText, htmlText);
        }

        public String extractText(Reader reader) {
            final ArrayList<String> list = new ArrayList<String>();

            HTMLEditorKit.ParserCallback parserCallback = new HTMLEditorKit.ParserCallback() {
                public void handleText(final char[] data, final int pos) {
                    list.add(new String(data));
                }

                public void handleStartTag(HTML.Tag tag, MutableAttributeSet attribute, int pos) {
                }

                public void handleEndTag(HTML.Tag t, final int pos) {
                }

                public void handleSimpleTag(HTML.Tag t, MutableAttributeSet a, final int pos) {
                    if (t.equals(HTML.Tag.BR)) {
                        list.add("\n");
                    }
                }

                public void handleComment(final char[] data, final int pos) {
                }

                public void handleError(final String errMsg, final int pos) {
                }
            };
            try {
                new ParserDelegator().parse(reader, parserCallback, true);
            } catch (IOException e) {
                e.printStackTrace();
            }
            String result = "";
            for (String s : list) {
                result += s;
            }
            return result;
        }


        @Override
        public void exportToClipboard(JComponent comp, Clipboard clip, int action) throws IllegalStateException {
            if (action == COPY) {
                clip.setContents(this.createTransferable(comp), null);
            }
        }

        @Override
        public int getSourceActions(JComponent c) {
            return COPY;
        }

    }

    class MyTransferable implements Transferable {


        private final String plainData;
        private final String htmlData;

        public MyTransferable(String plainData, String htmlData) {
            this.plainData = plainData;
            this.htmlData = htmlData;
        }

        public DataFlavor[] getTransferDataFlavors() {
            return supportedFlavors;
        }

        public boolean isDataFlavorSupported(DataFlavor flavor) {
            for (DataFlavor supportedFlavor : supportedFlavors) {
                if (supportedFlavor == flavor) {
                    return true;
                }
            }
            return false;
        }

        public Object getTransferData(DataFlavor flavor) throws UnsupportedFlavorException, IOException {
            if (flavor.equals(supportedFlavors[0])) {
                return htmlData;
            }
            if (flavor.equals(supportedFlavors[1])) {
                return plainData;
            }
            throw new UnsupportedFlavorException(flavor);
        }
    }


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

        TooltipTextFrame frame = new TooltipTextFrame("test", text.toString());
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


        TooltipTextFrame frame = new TooltipTextFrame("test", buf.toString());
        frame.setVisible(true);

    }
}
