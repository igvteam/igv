package org.broad.igv.feature.tribble;


import org.broad.tribble.readers.LineReader;
import org.broad.tribble.source.query.AsciiQuerySource;
import org.junit.Test;

import java.io.IOException;

/**
 * @author Jim Robinson
 * @date 10/12/11
 */
public class VCFWrapperCodecTest {

    @Test
    public void testGATKQuery() throws IOException {
        String path = "http://www.broadinstitute.org/igvdata/1KG/b37/GATK_bundle/current/b37/dbsnp_132.b37.vcf";

        AsciiQuerySource source = new AsciiQuerySource(path, path + ".idx");
        LineReader lineReader = source.query("17", 7569719, 7592864);
        String nextLine;
        while ((nextLine = lineReader.readLine()) != null) {
            System.out.println(nextLine);
        }
    }

}

