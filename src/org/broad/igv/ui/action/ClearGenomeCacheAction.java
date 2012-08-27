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

/*
* To change this template, choose Tools | Templates
* and open the template in the editor.
*/
package org.broad.igv.ui.action;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.log4j.Logger;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import java.awt.event.ActionEvent;

/**
 * @author eflakes
 */
public class ClearGenomeCacheAction extends MenuAction {

    static Logger logger = Logger.getLogger(ClearGenomeCacheAction.class);

    /**
     * Constructs ...
     *
     * @param label
     */
    public ClearGenomeCacheAction(String label) {
        super(label, null);
        setToolTipText(UIConstants.CLEAR_GENOME_CACHE_TOOLTIP);
    }

    /**
     * Method description
     *
     * @param evt
     */
    @Override
    public void actionPerformed(ActionEvent evt) {

        try {
            int option = (JOptionPane.showConfirmDialog(IGV.getMainFrame(),
                    "Clear the genome cache ?", "Clear the genome cache ?",
                    JOptionPane.YES_NO_OPTION));
            if (option == JOptionPane.YES_OPTION) {
                GenomeManager.getInstance().clearGenomeCache();
                UIUtilities.invokeOnEventThread(new Runnable() {
                    public void run() {
                        IGV.getInstance().rebuildGenomeDropdownList();
                    }
                });
            }


        } catch (Exception e) {
            JOptionPane.showMessageDialog(IGV.getMainFrame(),
                    "Error encontered while removing genomes: " + e.getMessage());
            logger.error("Error removing genomes from the user-defined genome list.", e);
        }
    }

}
