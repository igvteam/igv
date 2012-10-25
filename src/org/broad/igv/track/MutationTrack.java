package org.broad.igv.track;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.Mutation;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.TooltipTextFrame;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.NamedRunnable;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.Feature;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.net.URL;

/**
 * @author Jim Robinson
 * @date 11/8/11
 */
public class MutationTrack extends FeatureTrack {

    private static Logger log = Logger.getLogger(MutationTrack.class);

    public MutationTrack(ResourceLocator locator, String id, FeatureSource source) {
        super(locator, id, source);
        setSortable(true);
    }


    @Override
    public void overlay(RenderContext context, Rectangle rect) {
        if (!context.getChr().equals(Globals.CHR_ALL) || IGV.getInstance().getSession().getPreferenceAsBoolean(PreferenceManager.OVERLAY_MUTATIONS_WHOLE_GENOME)) {
            renderFeatures(context, rect);
        }
    }


    /**
     * Return a string for popup text.
     *
     * @param chr
     * @param position in genomic coordinates
     * @param y        - pixel position in panel coordinates (i.e. not track coordinates)
     * @return
     */
    @Override
    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {
        return super.getValueStringAt(chr, position, y, frame);    //To change body of overridden methods use File | Settings | File Templates.
    }

    @Override
    public boolean handleDataClick(TrackClickEvent te) {

        Feature f = getFeatureAtMousePosition(te);
        if (f != null && f instanceof Mutation) {

            final Mutation mut = (Mutation) f;
            final MouseEvent me = te.getMouseEvent();
            System.out.println("Submitting");
            LongRunningTask.submit(new NamedRunnable() {
                public String getName() {
                    return "Call OMA";
                }

                public void run() {

                    StringBuffer buf = new StringBuffer();
                    buf.append("<html>");

                    final String omaURL = mut.getOMAUrl();
                    if (omaURL != null) {
                        System.out.println("Running");
                        String omaText = getOMAText(mut, omaURL);
                        if (omaText != null) {
                            buf.append(omaText);
                            buf.append("<p>----------------------------------</p>");
                        }
                    }

                    buf.append(mut.getFullDescription());
                    buf.append("</html>");


                    final TooltipTextFrame tf = new TooltipTextFrame(MutationTrack.this.getName(), buf.toString());
                    Point p = me.getComponent().getLocationOnScreen();
                    tf.setLocation(Math.max(0, p.x + me.getX() - 150), Math.max(0, p.y + me.getY() - 150));

                    UIUtilities.invokeOnEventThread(new Runnable() {
                        public void run() {
                            tf.setVisible(true);
                        }
                    });

                }
            });


            return true;
        }
        return false;
    }


    private String getOMAText(Mutation mut, String url) {
        try {
            String omaWebService = url + "&frm=txt&fts=all";

            String result = HttpUtils.getInstance().getContentsAsString(new URL(omaWebService));

            BufferedReader br = new BufferedReader(new StringReader(result));
            String[] headers = br.readLine().split("\t");
            String[] values = br.readLine().split("\t");

            StringBuffer buf = new StringBuffer();
            buf.append("<p style=\"font-size:medium;\"><b>");
            buf.append(mut.getDescription());
            buf.append("</b><br><br>");

            buf.append("<table>");
            buf.append("<tr><td>Data from</td><td><a href=\"");
            buf.append(url);
            buf.append("\">Mutation Assessor</a></td><tr>");
            buf.append("<tr><td>Type</td><td>" + mut.getMutationType() + "</td></tr>");

            int n = Math.min(headers.length, values.length);
            for (int i = 0; i < n; i++) {
                final String header = headers[i].trim();
                final String value = values[i].trim();

                if (header.length() == 0 || value.length() == 0 || header.equals("Mutation") || header.equals("Type")) continue;

                buf.append("<tr>");
                buf.append("<td>");
                buf.append(header);
                buf.append("</td><td>");

                if (header.equals("MSA") || header.equals("PDB")) {
                    buf.append("<a href=\"" + value + "\">");
                    buf.append(header);
                    buf.append("</a>");
                } else if (header.equals("Uniprot")) {
                    buf.append("<a href=\"http://www.uniprot.org/uniprot/" + value + "\">");
                    buf.append(value);
                    buf.append("</a>");
                } else if (header.equals("Refseq")) {
                    buf.append("<a href=\"http://www.ncbi.nlm.nih.gov/sites/entrez?db=protein&cmd=search&term=" + value + "\">");
                    buf.append(value);
                    buf.append("</a>");
                } else if (header.toLowerCase().contains("impact")) {
                    if (value.toLowerCase().equals("high") || value.toLowerCase().equals("medium")) {
                        String color = value.toLowerCase().equals("high") ? "red" : "#0033FF";
                        buf.append("<div style=\"color:" + color + "\"><b>");
                        buf.append(value);
                        buf.append("</b></div>");
                    } else {
                        buf.append(value);
                    }
                } else {
                    buf.append(value);
                }
                buf.append("</td></tr>");
            }
            buf.append("</table></p>");

            return buf.toString();
        } catch (IOException e) {
            log.error("Error accessing OMA ", e);
            return null;
        }
    }

}
