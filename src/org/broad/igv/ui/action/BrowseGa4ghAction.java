/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.ui.action;

import org.apache.log4j.Logger;
import org.broad.igv.ga4gh.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.LongRunningTask;

import java.awt.event.ActionEvent;
import java.io.IOException;
import java.util.ArrayList;
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
                    final List<Ga4ghProvider> validProviders = new ArrayList<Ga4ghProvider>();

                    for (Ga4ghProvider provider : Ga4ghAPIHelper.providers) {
                        boolean valid = true;
                        for (Ga4ghDataset ds : provider.getDatasets()) {

                            List<Ga4ghReadset> readsets = Ga4ghAPIHelper.readsetSearch(provider, ds, 10);

                            if (readsets == null) {
                                log.error("No readsets returned for dataset: " + ds.getName());
                                valid = false;
                                // TODO -- for now, just exit, its confusing to bring up a partial list of providers
                                return;

                                //break;    // Something's wrong, probably authorization
                            }

                            ds.setReadsets(readsets);

                        }
                        if (valid) validProviders.add(provider);

                    }

                    if (validProviders.size() > 0) {
                        UIUtilities.invokeOnEventThread(new Runnable() {

                            @Override
                            public void run() {
                                Ga4ghLoadDialog dlg = (new Ga4ghLoadDialog(IGV.getMainFrame(), validProviders.toArray(new Ga4ghProvider[]{})));
                                dlg.setModal(true);
                                dlg.setVisible(true);
                                dlg.dispose();
                            }
                        });
                    }
                } catch (IOException e) {
                    log.error("Error opening Ga4gh dialog", e);
                }
            }
        });
    }

}
