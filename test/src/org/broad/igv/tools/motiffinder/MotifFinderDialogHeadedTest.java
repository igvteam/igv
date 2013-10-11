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

package org.broad.igv.tools.motiffinder;

import org.apache.commons.lang.StringUtils;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.ui.AbstractHeadedTest;
import org.broad.igv.ui.IGV;
import org.junit.Test;

import java.util.List;

import static junit.framework.Assert.assertEquals;

/**
 * @author jacob
 * @date 2013-Oct-11
 */
public class MotifFinderDialogHeadedTest extends AbstractHeadedTest{


    @Test
    public void testMultiLF() throws Exception{
        tstMultiPattern("\n");
    }

    @Test
    public void testMultiCRLF() throws Exception{
        tstMultiPattern("\r\n");
    }

    @Test
    public void testMultiCR() throws Exception{
        tstMultiPattern("\r");
    }


    private void tstMultiPattern(String sep) throws Exception{
        int initCount = IGV.getInstance().getFeatureTracks().size();

        String[] patterns = new String[]{"TATTAAT", "NCTC[GC]{3,8}", "NRTGC"};
        MotifFinderDialog dialog = new MotifFinderDialog(IGV.getMainFrame());
        String patternText = StringUtils.join(patterns, sep);

        dialog.setModal(false);
        dialog.setVisible(true);
        dialog.setPatternFieldText(patternText);

        dialog.okButtonActionPerformed(null);
        dialog.setVisible(false);

        MotifFinderPlugin.handleDialogResult(dialog);

        List<FeatureTrack> featureTrackList = IGV.getInstance().getFeatureTracks();
        assertEquals(initCount + patterns.length*2, featureTrackList.size());
    }
}
