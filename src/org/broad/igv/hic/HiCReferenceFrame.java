/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not
 * responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL), Version 2.1 which is
 * available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.hic;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ui.panel.ReferenceFrame;

/**
 * @author Jim Robinson
 * @date 5/3/12
 */
public class HiCReferenceFrame extends ReferenceFrame {

    Genome genome;

    public HiCReferenceFrame(String name, Genome genome) {
        super(name);
        this.genome = genome;
    }

    @Override
    protected Genome getGenome() {
        return genome;
    }
}
