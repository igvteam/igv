package org.igv.ui.action;


import org.igv.logging.*;
import org.igv.prefs.IGVPreferences;
import org.igv.prefs.PreferencesManager;
import org.igv.ui.IGV;
import org.igv.ui.UIConstants;
import org.igv.ui.util.MessageUtils;
import org.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.util.Arrays;

import static org.igv.prefs.Constants.*;

/**
 * @author jrobinso
 */
public class LoadFromDatabaseAction extends MenuAction {

    static Logger log = LogManager.getLogger(LoadFromDatabaseAction.class);
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

                String host = PreferencesManager.getPreferences().get(DB_HOST);
                if (host == null || host.trim().length() == 0) {
                    MessageUtils.showMessage("Please set database configuration in user preferences (View > Preferences)");
                    return null;
                }

                final IGVPreferences preferenceManager = PreferencesManager.getPreferences();
                String db = preferenceManager.get(DB_NAME);
                String port = preferenceManager.get(DB_PORT);

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
                loc1.setFormat(".seg");


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
