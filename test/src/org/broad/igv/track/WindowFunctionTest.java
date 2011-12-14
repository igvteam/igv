package org.broad.igv.track;

import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * Created by IntelliJ IDEA.
 * User: jacob
 * Date: 2011/01/05
 */
public class WindowFunctionTest {

    @Test
    public void testGetWindowName() throws Exception {

        WindowFunction wf = WindowFunction.percentile10;
        //System.out.println(wf.toString()); // "percentile10"
        //System.out.println(wf.getDisplayName());    // "10th Percentile"
        //System.out.println(wf.name());    // "percentile10"

        WindowFunction twf = WindowFunction.getWindowFunction("percentile10");

        assertEquals(wf, twf);
    }
}
