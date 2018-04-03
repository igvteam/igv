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

                            List<Ga4ghReadset> readsets = Ga4ghAPIHelper.searchReadGroupsets(provider, ds.getId(), 10);

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
