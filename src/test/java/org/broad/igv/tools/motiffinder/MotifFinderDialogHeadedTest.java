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
