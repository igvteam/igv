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
import org.broad.igv.util.ResourceLocator;

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
