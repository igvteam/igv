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

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

/**
 * @author jrobinso
 */
public class TileFactory {

    public static TDFTile createTile(byte[] buffer, int nSamples) throws IOException {

        ByteBuffer byteBuffer = ByteBuffer.wrap(buffer);
        byteBuffer.order(ByteOrder.LITTLE_ENDIAN);

        String typeString = StringUtils.readString(byteBuffer);
        TDFTile.Type type = TDFTile.Type.valueOf(typeString);

        switch (type) {
            case fixedStep:
                return new TDFFixedTile(byteBuffer, nSamples);
            case variableStep:
                return new TDFVaryTile(byteBuffer, nSamples);
            case bed:
            case bedWithName:
                return new TDFBedTile(byteBuffer, nSamples, type);
            default:
                throw new RuntimeException("Unknown tile type: " + type.toString());
        }
    }


}
