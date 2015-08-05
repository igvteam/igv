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

package org.broad.igv.gs.atm;

import java.util.List;

/**
* @author Jim Robinson
* @date Aug 3, 2011
*/
public class FileParameter {
    private String name;// 	Name of the parameter.
    private String required;// 	Boolean flag indicating whether or not the parameter is required.
    private String compositeFilename;//  	True means that the value for this parameter can contain multiple filenames.  The filenames are separated by the nameDelimeters below.
    private String nameDelimiters;// 	The string used to delimit multiple filenames.  This is only used if compositeFilename is true.
    private List<GSDataFormat> formats;// 	A set of GSDataFormats, describing the acceptable file formats for this FileParameter.  See below.

    public FileParameter(String fpName, String fpRequired, String fpCompositeFilename, String fpNameDelimiters,
                         List<GSDataFormat> dataFormats) {

        this.name = fpName;
        this.required = fpRequired;
        this.compositeFilename = fpCompositeFilename;
        this.nameDelimiters = fpNameDelimiters;
        this.formats = dataFormats;
    }

    public void print() {
        System.out.print(name + "  ");
        for(GSDataFormat df : formats) {
            System.out.print(" " + df + ",");
        }
        System.out.println();

    }

    public String getName() {
        return name;
    }

    public String getRequired() {
        return required;
    }

    public String getCompositeFilename() {
        return compositeFilename;
    }

    public String getNameDelimiters() {
        return nameDelimiters;
    }

    public List<GSDataFormat> getFormats() {
        return formats;
    }
}
