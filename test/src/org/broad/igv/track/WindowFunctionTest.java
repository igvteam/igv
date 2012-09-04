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
