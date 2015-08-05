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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.tdf;

import org.broad.igv.util.StringUtils;

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
