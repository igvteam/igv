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

package org.broad.igv.session;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.igv.util.Utilities;
import org.junit.Before;
import org.junit.Test;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

import java.io.File;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012-Aug-07
 */
public class IGVSessionReaderTest extends AbstractHeadlessTest {

    IGVSessionReader sessionReader;

    @Before
    public void setUp() throws Exception {
        super.setUp();
        sessionReader = new IGVSessionReader(null);

    }

    /**
     * Read a session file which uses relative paths, check that the paths
     * are converted to absolute and loaded properly
     *
     * @throws Exception
     */
    @Test
    public void testReadRelativePaths() throws Exception {
        String sessionPath = TestUtils.DATA_DIR + "sessions/testBedsRelPath.xml";
        Session session = new Session(sessionPath);

        InputStream inputStream = ParsingUtils.openInputStreamGZ(new ResourceLocator(sessionPath));
        Document document = Utilities.createDOMDocumentFromXmlStream(inputStream);
        NodeList elements = document.getElementsByTagName("Resource");

        sessionReader.dataFiles = new ArrayList<ResourceLocator>();
        for (int el = 0; el < elements.getLength(); el++) {
            Element element = (Element) elements.item(el);
            sessionReader.processResource(session, element, new HashMap());
        }
        assertEquals(elements.getLength(), sessionReader.dataFiles.size());

        TrackLoader loader = new TrackLoader();
        for (ResourceLocator locator : sessionReader.dataFiles) {
            File file = new File(locator.getPath());
            assertTrue(file.isAbsolute());
            assertTrue(file.canRead());

            File testDataDir = new File(TestUtils.DATA_DIR);
            assertTrue(file.getAbsolutePath().contains(testDataDir.getAbsolutePath()));

            List<Track> testTrack = loader.load(locator, genome);
            assertEquals(1, testTrack.size());
        }
    }
}
