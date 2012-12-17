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
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.broad.igv.util.Utilities;
import org.junit.Test;
import org.w3c.dom.Document;

import static org.junit.Assert.assertEquals;

/**
 * User: jacob
 * Date: 2012-Dec-17
 */
public class RecursiveAttributesTest extends AbstractHeadlessTest {
    @Test
    public void testReadElement() throws Exception {
        String sessionPath = TestUtils.DATA_DIR + "sessions/BRCA_loh.xml";
        Document document = Utilities.createDOMDocumentFromXmlStream(ParsingUtils.openInputStreamGZ(new ResourceLocator(sessionPath)));

        RecursiveAttributes allAttributes = RecursiveAttributes.readElement(document.getDocumentElement());

        assertEquals("Session", allAttributes.getName());

        assertEquals(3, allAttributes.getAttributes().size());
        assertEquals("hg18", allAttributes.get("genome"));

        assertEquals(4, allAttributes.getChildren().size());

        int panelCount = 0;
        for(RecursiveAttributes childAttributes: allAttributes.getChildren()){
            if(childAttributes.getName().equalsIgnoreCase("Panel")){
                panelCount++;
                if(childAttributes.get("name").equals("FeaturePanel")){
                    assertEquals(2, childAttributes.getChildren().size());
                }
            }
        }
        assertEquals(2, panelCount);
    }
}
