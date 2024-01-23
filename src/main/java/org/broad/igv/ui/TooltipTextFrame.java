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

import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.util.BrowserLauncher;
import org.broad.igv.util.HttpUtils;

import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import javax.swing.*;
import javax.swing.event.HyperlinkEvent;

/**
 * @author Jim Robinson
 * @date 11/7/11
 */
public class TooltipTextFrame extends JFrame {

    private static Logger log = LogManager.getLogger(TooltipTextFrame.class);

    public TooltipTextFrame(String title, String text) throws HeadlessException {

        setTitle(title);

        setAlwaysOnTop(true);

        setLayout(new BorderLayout());

        // Translate line breaks to divs.  This is a workaround to allow partial selection of text and maintain
        // line feeds.  The HTMLDocument does not preserve line feeds in text output for selectedText with <br> tags,
        // but does with <div>
        text = text
                .replace("<br>", "<div>")
                .replace("<br/>", "<div>");

        JEditorPane pane = new JEditorPane("text/html", text);
        pane.setEditable(false);
        pane.setTransferHandler(new TTTransferHandler());

        Dimension d = pane.getPreferredSize();
        int w = (int) (1.2 * d.width);
        int h = (int) (1.25 * d.height);

        h = h > 600 ? 600 : (Math.max(h, 100));
        w = w > 800 ? 800 : (Math.max(w, 100));
        setSize(w, h);

        JScrollPane scrollPane = new JScrollPane(pane);
        pane.addHyperlinkListener(e -> {
            try {
                if (e.getEventType() == HyperlinkEvent.EventType.ACTIVATED)
                    BrowserLauncher.openURL(e.getURL().toExternalForm());
            } catch (IOException e1) {
                log.error("Error opening hyperlink", e1);
            }
        });

        getContentPane().add(scrollPane, BorderLayout.CENTER);
    }


    /**
     * TransferHandler for copy to clipboard -- fetches selected text from the pane, if no text is selected
     * fetches all text.
     */
    static class TTTransferHandler extends TransferHandler {

        protected Transferable createTransferable(JComponent c) {
            final JEditorPane pane = (JEditorPane) c;
            String selectedText = pane.getSelectedText();
            if (selectedText == null) {
                pane.selectAll();
                selectedText = pane.getSelectedText();
            }
            return new TTTransferable(selectedText);
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

    static class TTTransferable implements Transferable {

        private static final DataFlavor[] supportedFlavors;

        static {
            try {
                supportedFlavors = new DataFlavor[]{
                        new DataFlavor("text/plain;class=java.lang.String")
                };
            } catch (ClassNotFoundException e) {
                throw new ExceptionInInitializerError(e);
            }
        }

        private final String text;

        public TTTransferable(String text) {
            this.text = text;
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

        public Object getTransferData(DataFlavor flavor) throws UnsupportedFlavorException {
            if (flavor.equals(supportedFlavors[0])) {
                return text;
            } else {
                throw new UnsupportedFlavorException(flavor);
            }
        }
    }


    public static void main(String[] args) throws IOException {

        test1();
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
        String result = HttpUtils.getInstance().getContentsAsString(HttpUtils.createURL(url));

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
