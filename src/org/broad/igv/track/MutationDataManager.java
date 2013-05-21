/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.track;

import org.broad.igv.feature.Mutation;
import org.broad.igv.feature.genome.Genome;

import java.io.IOException;

/**
 * @author jrobinso
 *         Date: 4/9/13
 *         Time: 8:49 AM
 */
public class MutationDataManager extends MultitrackDataManager<Mutation>{

    public MutationDataManager(String path, Genome genome) throws IOException {
        super(path, genome);
    }

    @Override
    protected String getTrackKey(Mutation feat) {
        return feat.getSampleId();
    }


}
