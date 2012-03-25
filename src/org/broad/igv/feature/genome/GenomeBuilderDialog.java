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

package org.broad.igv.feature.genome;

import org.broad.igv.Globals;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.OkCancelDialog;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.File;


/**
 * @author eflakes
 */
public class GenomeBuilderDialog extends OkCancelDialog {

    private GenomeBuilderPane builderPane;
    private File genomeArchiveFile = null;

    public GenomeBuilderDialog(IGV igv, boolean modal) {

        super(igv.getMainFrame(), modal);
        builderPane = new GenomeBuilderPane(igv);
        setTitle("Import Genome");

        setSize(800, 500);
        setResizable(false);

        JPanel contentPane = getDialogPanel();
        contentPane.setLayout(new BorderLayout());
        contentPane.add(builderPane, BorderLayout.CENTER);
        setLocationRelativeTo(igv.getMainFrame());
        setOkButtonText(" Save ");
    }

    public String getCytobandFileName() {
        return builderPane.getCytobandFileName();
    }

    public String getFastaFileName() {
        return builderPane.getFastaFileName();
    }

    public String getChrAliasFileName() {
        return builderPane.getChrAliasFileName();
    }

    public String getGenomeId() {
        return builderPane.getGenomeId();
    }

    public String getGenomeDisplayName() {
        return builderPane.getGenomeDisplayName();
    }

    public String getRefFlatFileName() {
        return builderPane.getRefFlatFileName();
    }

    public String getGenomeArchiveLocation() {
        return builderPane.getGenomeArchiveLocation();
    }

    public String getArchiveFileName() {
        return builderPane.getArchiveFileName();
    }

    public String getSequenceLocation() {
        return getFastaFileName();
    }

    public String getSequenceLocationOverride() {
        String seqLocation = builderPane.getSequenceURL();
        if (seqLocation == null && seqLocation.trim().length() == 0) {
            return null;
        }
        return seqLocation;

    }

    @Override
    public boolean cancelButtonClicked(ActionEvent event) {
        return true;
    }

    @Override
    public boolean okButtonClicked(ActionEvent event) {

        boolean isOk = builderPane.validateSelection();

        // Passed validation now get genome location and check it
        if (isOk) {
            genomeArchiveFile = builderPane.showGenomeArchiveDirectoryChooser();
            if (genomeArchiveFile == null) {
                isOk = false;
            }
        }
        return isOk;
    }
}