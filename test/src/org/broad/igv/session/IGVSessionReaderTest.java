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

package org.broad.igv.session;

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.renderer.Renderer;
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

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.adapters.XmlJavaTypeAdapter;
import java.awt.*;
import java.io.File;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import static junit.framework.Assert.*;

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
            sessionReader.processResource(session, element, new HashMap(), sessionPath, null);
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

    @Test
    public void testReadWriteColorString() throws Exception{
        SessionDataHolder inVal = new SessionDataHolder();
        Color inCol = new Color(0, 100, 200);
        inVal.color = inCol;
        SessionDataHolder outVal = TestUtils.marshallUnmarshall(inVal);
        assertEquals(inCol, outVal.color);
        assertNull(outVal.renderer);
    }

    @Test
    public void testReadWriteRenderer() throws Exception{
        SessionDataHolder inVal = new SessionDataHolder();
        inVal.renderer = (Renderer) RendererFactory.defaultRendererClass.newInstance();

        assert inVal.renderer != null;

        SessionDataHolder outVal = TestUtils.marshallUnmarshall(inVal);
        assertEquals(inVal.renderer.getClass().getName(), outVal.renderer.getClass().getName());
    }

    private static class SessionDataHolder{
        @XmlJavaTypeAdapter(SessionXmlAdapters.Color.class)
        @XmlAttribute public Color color;

        @XmlJavaTypeAdapter(SessionXmlAdapters.Renderer.class)
        @XmlAttribute public Renderer renderer;
    }


}
