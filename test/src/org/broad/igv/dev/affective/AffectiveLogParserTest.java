package org.broad.igv.dev.affective;

import org.junit.Test;

import java.io.*;
import java.text.ParseException;

import static org.junit.Assert.assertEquals;

/**
 * @author Jim Robinson
 * @date 1/22/12
 */
public class AffectiveLogParserTest {

    String testFile = "test/data/affective/test.affectiva.csv";


    @Test
    public void testParseHeader() throws IOException, ParseException {

        // Start Time: 2011-04-06 08:57:35 Offset:-04
        // Z-axis | Y-axis | X-axis | Battery | °Celsius | EDA(uS)
        String expectedStartDate = "2011-04-06";
        int expectedStartTime = (8 * 60 + 57) * 60 + 35;
        String[] expectedTrackNames = {"Z-axis", "Y-axis", "X-axis", "Battery", "°Celsius", "EDA(uS)"};
        int expectedSamplingRate = 8;
        String expectedUUID = "ABC441000XX";

        AffectiveLogParser parser = new AffectiveLogParser(testFile, null);
        BufferedReader reader = new BufferedReader(new FileReader(testFile));
        AffectiveLogParser.Header header = parser.parseHeader(reader);
        reader.close();

        assertEquals("Start date", expectedStartDate, header.getStartDate());
        assertEquals("Start time", expectedStartTime, header.getStartTime());
        assertEquals("UUID", expectedUUID, header.getUuid());
        assertEquals("Sampling rate", expectedSamplingRate, header.getSamplingRate());
        assertEquals(expectedTrackNames.length, header.getTrackNames().length);
        for (int i = 0; i < expectedTrackNames.length; i++) {
            assertEquals("Track name", expectedTrackNames[i], header.getTrackNames()[i]);
        }

    }

}
