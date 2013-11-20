package org.broad.igv.ui.action;

import org.apache.log4j.Logger;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.session.SessionXmlAdapters;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.collections.CI;
import org.broad.igv.util.encode.EncodeFileBrowser;
import org.broad.igv.util.encode.EncodeFileRecord;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 *         Date: 11/2/13
 *         Time: 6:39 PM
 */
public class BrowseEncodeAction extends MenuAction {

    private static Logger log = Logger.getLogger(BrowseEncodeAction.class);

    private static Map<String, Color> colors;

    static {
        colors = new HashMap<String, Color>();
        colors.put("H3K27AC", new Color(200, 0, 0));
        colors.put("H3K27ME3", new Color(200, 0, 0));
        colors.put("H3K36ME3", new Color(0, 0, 150));
        colors.put("H3K4ME1", new Color(0, 150, 0));
        colors.put("H3K4ME2", new Color(0, 150, 0));
        colors.put("H3K4ME3", new Color(0, 150, 0));
        colors.put("H3K9AC", new Color(100, 0, 0));
        colors.put("H3K9ME1", new Color(100, 0, 0));
    }


    IGV igv;

    public BrowseEncodeAction(String label, int mnemonic, IGV igv) {
        super(label, null, mnemonic);
        this.igv = igv;
    }


    @Override
    public void actionPerformed(ActionEvent event) {

        try {
            Genome genome = igv.getGenomeManager().getCurrentGenome();
            EncodeFileBrowser browser = EncodeFileBrowser.getInstance(genome.getId());

            if (browser == null) {
                MessageUtils.showMessage("Encode data is not available for " + genome.getDisplayName() + " through IGV.");
                return;
            }

            browser.setVisible(true);
            if (browser.isCanceled()) return;

            java.util.List<EncodeFileRecord> records = browser.getSelectedRecords();
            if (records.size() > 0) {
                List<ResourceLocator> locators = new ArrayList<ResourceLocator>(records.size());
                for (EncodeFileRecord record : records) {
                    ResourceLocator rl = new ResourceLocator(record.getPath());
                    rl.setName(record.getTrackName());

                    final String antibody = record.getAttributeValue("antibody");
                    if (antibody != null) {
                        rl.setColor(colors.get(antibody.toUpperCase()));
                    }

                    locators.add(rl);
                }
                igv.loadTracks(locators);
            }


        } catch (IOException e) {
            log.error("Error opening Encode browser", e);
        }
    }
}
