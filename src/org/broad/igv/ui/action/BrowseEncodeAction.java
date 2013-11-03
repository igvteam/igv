package org.broad.igv.ui.action;

import org.apache.log4j.Logger;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.encode.EncodeFileBrowser;
import org.broad.igv.util.encode.EncodeFileRecord;

import java.awt.event.ActionEvent;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * @author jrobinso
 *         Date: 11/2/13
 *         Time: 6:39 PM
 */
public class BrowseEncodeAction extends MenuAction {

    private static Logger log = Logger.getLogger(BrowseEncodeAction.class);

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

            if(browser == null) {
                MessageUtils.showMessage("Encode data is not available for " + genome.getDisplayName() + " through IGV.");
            }

            browser.setVisible(true);
            if(browser.isCanceled()) return;

           java.util.List<EncodeFileRecord> records = browser.getSelectedRecords();
            if(records.size() > 0) {
                List<ResourceLocator> locators = new ArrayList<ResourceLocator>(records.size());
                for(EncodeFileRecord record : records) {
                    locators.add(new ResourceLocator(record.getPath()));
                }
                igv.loadTracks(locators);
            }


        } catch (IOException e) {
            log.error("Error opening Encode browser", e);
        }
    }
}
