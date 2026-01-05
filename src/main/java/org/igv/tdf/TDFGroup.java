package org.igv.tdf;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.Map;

/**
 * @author jrobinso
 */
public class TDFGroup extends TDFEntity {

    public static final String USE_PERCENTILE_AUTOSCALING = "userPercentileAutoscaling";


    TDFGroup(String name) {
        super(name);
    }

    TDFGroup(String name, Map<String, String> attributes) {
        super(name, attributes);
    }

    public TDFGroup(String name, ByteBuffer byteBuffer) throws IOException {
        super(name);
        fill(byteBuffer);
    }

    void write(BufferedByteWriter dos) throws IOException {

        writeAttributes(dos);
    }

    void fill(ByteBuffer byteBuffer) throws IOException {
        readAttributes(byteBuffer);
    }
}
