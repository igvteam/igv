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

package org.broad.igv.ui.util;

/**
 * @author Jim Robinson
 * @date 1/25/12
 */
public class ImageFileTypes {

    /**
     * Snapshot types
     */
    public static enum Type {

        NULL("", ""),
        EPS(".eps", "Encapsulated Postscript Files (*.eps)"),
        PDF(".pdf", "Portable Document FormatFles (*.pdf)"),
        SVG(".svg", "Scalable Vector Graphics Files (*.svg)"),
        PNG(".png", "Portable Network Graphics Files (*.png)"),
        JPEG(".jpeg", "Joint Photographic Experts Group Files (*.jpeg)");
        private String fileExtension;
        private String fileDescription;

        Type(String extension, String description) {
            fileExtension = extension;
            fileDescription = description;
        }

        public String getExtension() {
            return fileExtension;
        }

        public String getDescription() {
            return fileDescription;
        }
    }

    public static Type getImageFileType(String fileExtension) {

        String extension = fileExtension.toLowerCase();
        Type type = Type.NULL;

        if(".jpg".equalsIgnoreCase(fileExtension)){
            type = Type.JPEG;
        }

        for(Type iterType: Type.values()){
            if(type != Type.NULL) break;
            if(iterType.getExtension().equalsIgnoreCase(extension)){
                type = iterType;
            }
        }

        return type;
    }

}
