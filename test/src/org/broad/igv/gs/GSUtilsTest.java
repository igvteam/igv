package org.broad.igv.gs;

import org.junit.Test;

import java.net.URL;

import static junit.framework.Assert.assertEquals;

/**
 * @author jrobinso
 *         Date: 7/22/13
 *         Time: 1:00 PM
 */
public class GSUtilsTest {

    @Test
    public void testParseDataFormatString() throws Exception {

        String dataFormatString = "http://www.genomespace.org/datamanager/dataformat/gct/0.0.0";

        String ext = GSUtils.parseDataFormatString(dataFormatString);

        assertEquals(".gct", ext);




    }
}
