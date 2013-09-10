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

package org.broad.igv.plugin.mongocollab;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.File;

import static org.junit.Assert.assertEquals;

/**
 * @author jacob
 * @date 2013-Sep-10
 */
public class MongoCollabPluginTest extends AbstractHeadlessTest {

    @Test
    public void testReadSpec() throws Exception{
        String filePath = TestUtils.DATA_DIR + "testMongo.txt";
        MongoCollabPlugin.Locator locator = new MongoCollabPlugin.Locator(new File(filePath));
        assertEquals("localhost", locator.host);
        assertEquals(27017, locator.port);
    }
}
