/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.action;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.NamedRunnable;
import org.broad.igv.util.ResourceLocator;
import sun.util.resources.CurrencyNames_vi_VN;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.util.Arrays;

/**
 * @author jrobinso
 */
public class LoadFromDatabaseAction extends MenuAction {

    static Logger log = Logger.getLogger(LoadFromDatabaseAction.class);
    IGV mainFrame;


    public LoadFromDatabaseAction(String label, int mnemonic, IGV mainFrame) {
        super(label, null, mnemonic);
        this.mainFrame = mainFrame;
        setToolTipText(UIConstants.LOAD_SERVER_DATA_TOOLTIP);
    }


    @Override
    public void actionPerformed(ActionEvent evt) {


        SwingWorker worker = new SwingWorker() {

            @Override
            protected Object doInBackground() throws Exception {

                String host = PreferenceManager.getInstance().get(PreferenceManager.DB_HOST);
                if (host == null || host.trim().length() == 0) {
                    MessageUtils.showMessage("Please set database configuration in user preferences (View > Preferences)");
                    return null;
                }

                final PreferenceManager preferenceManager = PreferenceManager.getInstance();
                String db = preferenceManager.get(PreferenceManager.DB_NAME);
                String port = preferenceManager.get(PreferenceManager.DB_PORT);

                String url = "jdbc:mysql://" + host;
                if (!port.equals("-1")) {
                    url += ":" + port;
                }
                url += "/" + db;

                String table2 = "SAMPLE_INFO";
                ResourceLocator loc2 = new ResourceLocator(url, table2);
                loc2.setDescription("SELECT * FROM " + table2);
//

                String table1 = "CNV";
                ResourceLocator loc1 = new ResourceLocator(url, table1);
                // TODO -- get these mappings from a config table
                String query = "SELECT  Sample Sample, `Probe Median` Value, " +
                        "Chromosome chr,  Start start, Stop end, " +
                        "CONCAT('<br>Event: ', Event,'<br>% CNV Overlap = ', `% of CNV Overlap`) description " +
                        " FROM CNV " +
                        " WHERE Event = 'CN Gain' OR Event = 'CN Loss' OR Event = 'High Copy Gain'";
//                         + "INNER JOIN SAMPLE_INFO ON SAMPLE_INFO.SAMPLE = CNV.SAMPLE " +
//                        "WHERE SAMPLE_INFO.SUBTYPE like 'Classical'";

//                String query = "select * from cnv";
                loc1.setDescription(query);
                loc1.setType(".seg");


                mainFrame.loadTracks(Arrays.asList(loc1, loc2));

                return null;
            }

            @Override
            protected void done() {
                mainFrame.showLoadedTrackCount();
            }
        };

        worker.execute();


    }


}
