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
import org.broad.igv.feature.genome.GenomeListItem;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.IGV;
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
    }

    /**
     * Method description
     *
     * @param e
     */
    @Override
    public void actionPerformed(ActionEvent e) {

        try {

            Collection<GenomeListItem> genomeItemList =
                    GenomeManager.getInstance().getUserDefinedGenomeArchiveList();

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
                    GenomeManager.getInstance().updateImportedGenomePropertyFile();
                    IGV.getInstance().rebuildGenomeDropdownList();
                }
            }
        } catch (Exception ex) {
            JOptionPane.showMessageDialog(IGV.getMainFrame(),
                    "Error encontered while removing genomes: " + ex.getMessage());
            logger.error("Error removing genomes from the imported genome list.", ex);
        }
    }

}
