package org.broad.igv.track;

import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * Created by IntelliJ IDEA.
 * User: jacob
 * Date: 2011/12/14
 */
public class WindowFunctionTest {

    @Test
    public void testGetWindowName() throws Exception {

        WindowFunction wf = WindowFunction.percentile10;

        WindowFunction twf = WindowFunction.getWindowFunction("10th Percentile");

        assertEquals(wf, twf);
    }
}
