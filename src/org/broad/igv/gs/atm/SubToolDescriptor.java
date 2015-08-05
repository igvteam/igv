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
public class SubToolDescriptor {
    String name;//	Name of the web tool.
    String id;//	 	LSID which uniquely identifies the web tool.
    String version;//	 	The version of the web tool.
    String author;//	 	The author of the web tool.
    String description;//	 	Description of the web tool.
    String help;//	 	URI for help about the web tool.
    List<FileParameter> fileParameters;//		A list of FileParameters.  See below.
    String urlModifier; // 	The URL fragment that will be appended to the parent WebTools dmServer.

    public SubToolDescriptor(String name, String id, String version, String author, String description, String help,
                             String urlModifier, List<FileParameter> fileParameters) {
        this.name = name;
        this.id = id;
        this.version = version;
        this.author = author;
        this.description = description;
        this.help = help;
        this.fileParameters = fileParameters;
        this.urlModifier = urlModifier;
    }

    public void print() {
        System.out.println();
        System.out.println(name);
        for (FileParameter fp : fileParameters) {
            fp.print();
        }

    }
}
