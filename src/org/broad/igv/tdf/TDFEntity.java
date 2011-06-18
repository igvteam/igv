/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
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
