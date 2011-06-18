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
                IGV.getInstance().getGenomeManager().clearGenomeCache();
                UIUtilities.invokeOnEventThread(new Runnable() {
                    public void run() {
                        IGV.getInstance().rebuildGenomeDropdownList(null);
                    }
                });
            }


        }
        catch (Exception e) {
            JOptionPane.showMessageDialog(IGV.getMainFrame(),
                    "Error encontered while removing genomes: " + e.getMessage());
            logger.error("Error removing genomes from the user-defined genome list.", e);
        }
    }

}
