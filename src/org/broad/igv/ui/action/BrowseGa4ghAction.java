package org.broad.igv.ui.action;

import org.apache.log4j.Logger;
import org.broad.igv.ga4gh.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.LongRunningTask;

import java.awt.event.ActionEvent;
import java.io.IOException;
import java.util.List;

/**
 * Created by jrobinso on 8/16/14.
 */
public class BrowseGa4ghAction extends MenuAction {

    private static Logger log = Logger.getLogger(BrowseGa4ghAction.class);


    IGV igv;

    public BrowseGa4ghAction(String label, int mnemonic, IGV igv) {
        super(label, null, mnemonic);
        this.igv = igv;
    }


    @Override
    public void actionPerformed(ActionEvent event) {

        LongRunningTask.submit(new Runnable() {
            public void run() {
                try {
                    for (Ga4ghProvider provider : Ga4ghAPIHelper.providers) {

                        for (Ga4ghDataset ds : provider.getDatasets()) {

                            List<Ga4ghReadset> readsets = Ga4ghAPIHelper.readsetSearch(provider, ds, 10);

                            ds.setReadsets(readsets);

                        }
                    }

                    UIUtilities.invokeOnEventThread(new Runnable() {

                        @Override
                        public void run() {
                            Ga4ghLoadDialog dlg = (new Ga4ghLoadDialog(IGV.getMainFrame(), Ga4ghAPIHelper.providers));
                            dlg.setModal(true);
                            dlg.setVisible(true);
                            dlg.dispose();
                        }

                    });
                } catch (IOException e) {
                    log.error("Error opening Ga4gh dialog", e);
                }
            }
        });
    }

}
