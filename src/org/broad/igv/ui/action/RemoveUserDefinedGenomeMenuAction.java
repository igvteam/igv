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
import org.broad.igv.feature.genome.GenomeManager.GenomeListItem;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.UserDefinedGenomeCheckList;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.util.*;

/**
 * @author eflakes
 */
public class RemoveUserDefinedGenomeMenuAction extends MenuAction {

    static Logger logger = Logger.getLogger(RemoveUserDefinedGenomeMenuAction.class);

    /**
     * Constructs ...
     *
     * @param label
     * @param mnemonic
     */
    public RemoveUserDefinedGenomeMenuAction(String label, int mnemonic) {
        super(label, null, mnemonic);
        setToolTipText(UIConstants.REMOVE_USER_DEFINE_GENOME_TOOLTIP);
    }

    /**
     * Method description
     *
     * @param e
     */
    @Override
    public void actionPerformed(ActionEvent e) {

        try {

            List<GenomeListItem> genomeItemList =
                    IGV.getInstance().getGenomeManager().getUserDefinedGenomeArchiveList();

            if (genomeItemList.isEmpty()) {
                JOptionPane.showMessageDialog(IGV.getMainFrame(),
                        "There are no imported genomes to remove.");
            } else {

                GenomeListItem currentlySelectedDropdownGenome = IGV.getInstance().getGenomeSelectedInDropdown();
                List<String> genomeNames = new ArrayList();
                UserDefinedGenomeCheckList checkList = new UserDefinedGenomeCheckList(false);
                for (GenomeListItem genomeListItem : genomeItemList) {
                    genomeNames.add(genomeListItem.getDisplayableName());
                }
                checkList.addItems(genomeNames);
                checkList.sort();
                int status = JOptionPane.showConfirmDialog(IGV.getMainFrame(), checkList,
                        "Imported Genomes to Remove", JOptionPane.OK_CANCEL_OPTION,
                        JOptionPane.PLAIN_MESSAGE, null);

                if ((status == JOptionPane.CANCEL_OPTION) || (status == JOptionPane.CLOSED_OPTION)) {
                    return;
                }

                boolean removed = false;
                HashSet<String> selectedGenomes = checkList.getSelectedGenomes();
                Iterator iterator = genomeItemList.iterator();
                while (iterator.hasNext()) {
                    String genomeName = ((GenomeListItem) iterator.next()).getDisplayableName();

                    // Skip genome if not selected for removal
                    if (!selectedGenomes.contains(genomeName)) {
                        continue;
                    }

                    if (currentlySelectedDropdownGenome.getDisplayableName().equalsIgnoreCase(genomeName)) {
                        JOptionPane.showMessageDialog(
                                IGV.getMainFrame(),
                                "<html>Genome [" + genomeName + "] is currently in use and cannot be removed." +
                                        "<br>Please select another genome to view before trying to remove it.</html>");
                        continue;
                    }

                    if (selectedGenomes.contains(genomeName)) {
                        iterator.remove();
                        removed = true;
                    }
                }

                if (removed) {
                    IGV.getInstance().getGenomeManager().updateImportedGenomePropertyFile();
                    IGV.getInstance().rebuildGenomeDropdownList(null);
                }
            }
        }
        catch (Exception ex) {
            JOptionPane.showMessageDialog(IGV.getMainFrame(),
                    "Error encontered while removing genomes: " + ex.getMessage());
            logger.error("Error removing genomes from the imported genome list.", ex);
        }
    }

}
