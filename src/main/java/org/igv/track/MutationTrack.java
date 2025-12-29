package org.igv.track;

import htsjdk.tribble.Feature;
import org.igv.logging.*;
import org.igv.Globals;
import org.igv.feature.Mutation;
import org.igv.prefs.Constants;
import org.igv.ui.IGV;
import org.igv.ui.TooltipTextFrame;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.ui.util.UIUtilities;
import org.igv.util.HttpUtils;
import org.igv.util.LongRunningTask;
import org.igv.util.NamedRunnable;
import org.igv.util.ResourceLocator;

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

    private static Logger log = LogManager.getLogger(MutationTrack.class);

    public MutationTrack(ResourceLocator locator, String id, FeatureSource source) {
        super(locator, id, source);
        setSortable(true);
    }

    public MutationTrack() {
    }

    @Override
    public void overlay(RenderContext context, Rectangle rect) {
        if (!context.getChr().equals(Globals.CHR_ALL) ||
                IGV.getInstance().getSession().getPreferenceAsBoolean(Constants.OVERLAY_MUTATIONS_WHOLE_GENOME)) {
            renderFeatures(context, rect);
        }
    }


    /**
     * Return a string for popup text.
     *
     * @param chr
     * @param position in genomic coordinates
     * @param mouseX
     * @return
     */
    @Override
    public String getValueStringAt(String chr, double position, int mouseX, int mouseY, ReferenceFrame frame) {
        return super.getValueStringAt(chr, position, mouseX, mouseY, frame);    //To change body of overridden methods use File | Settings | File Templates.
    }

    @Override
    public boolean handleDataClick(TrackClickEvent te) {
        super.openTooltipWindow(te);
        return true;
    }


    private String getOMAText(Mutation mut, String url) {
        try {
            String omaWebService = url + "&frm=txt&fts=all";

            String result = HttpUtils.getInstance().getContentsAsString(HttpUtils.createURL(omaWebService));

            BufferedReader br = new BufferedReader(new StringReader(result));
            String[] headers = br.readLine().split("\t");
            String[] values = br.readLine().split("\t");

            StringBuffer buf = new StringBuffer();

            buf.append("<table>");
            buf.append("<tr><td>Data from</td><td><a href=\"");
            buf.append(url);
            buf.append("\">Mutation Assessor</a></td><tr>");
            buf.append("<tr><td>Type</td><td>" + mut.getMutationType() + "</td></tr>");

            int n = Math.min(headers.length, values.length);
            for (int i = 0; i < n; i++) {
                final String header = headers[i].trim();
                final String value = values[i].trim();

               if (header.length() == 0 || value.length() == 0) continue;

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
