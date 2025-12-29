package org.broad.igv.tools.motiffinder;

import org.apache.commons.lang3.StringUtils;
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
        MotifFinderDialog dialog = new MotifFinderDialog(IGV.getInstance().getMainFrame());
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
