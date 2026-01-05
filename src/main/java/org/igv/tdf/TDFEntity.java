package org.igv.tdf;

import org.igv.util.StringUtils;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * @author jrobinso
 */
public class TDFEntity {

    Map<String, String> attributes;

    private String name;

    public TDFEntity(String name) {
        this.name = name;
        this.attributes = new HashMap();
    }

    public TDFEntity(String name, Map<String, String> attributes) {
        this.name = name;
        this.attributes = attributes;
    }

    public Set<String> getAttributeNames() {
        return attributes.keySet();
    }

    public String getAttribute(String name) {
        return attributes.get(name);
    }

    public void setAttribute(String name, String value) {
        attributes.put(name, value);
    }

    void writeAttributes(BufferedByteWriter dos) throws IOException {
        dos.putInt(attributes.size());
        for (Map.Entry<String, String> entry : attributes.entrySet()) {
            writeString(dos, entry.getKey());
            writeString(dos, entry.getValue());
        }
    }

    void readAttributes(ByteBuffer byteBuffer) throws IOException {
        int nAttributes = byteBuffer.getInt();
        while (nAttributes-- > 0) {
            String key = StringUtils.readString(byteBuffer);
            String value = StringUtils.readString(byteBuffer);
            setAttribute(key, value);
        }
    }

    /**
     * Write a string as a null terminated character sequence.
     * <p/>
     * IGV requires all strings to be ascii,  so single byte storaged is enough
     *
     * @param dos
     * @param s
     * @throws java.io.IOException
     */
    public void writeString(BufferedByteWriter dos, String s) throws IOException {
        dos.putNullTerminatedString(s);
    }

    /**
     * @return the name
     */
    public String getName() {
        return name;
    }
}
