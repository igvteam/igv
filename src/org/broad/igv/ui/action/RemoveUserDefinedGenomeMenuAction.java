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
