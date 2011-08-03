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

package org.broad.igv.gs;

import java.util.List;
import java.util.Set;

/**
 * @author Jim Robinson
 * @date Aug 2, 2011
 */
public class WebToolDescriptor {
    String name;
    String id;
    String version;
    String author;
    String description;
    String help;
    List<FileParameter> fileParameters;
    List<SubToolDescriptor> subTools;
    String baseUrl;

    public WebToolDescriptor(String name, String id, String version, String author, String description, String help,
                             String baseUrl, List<FileParameter> fileParameters, List<SubToolDescriptor> subTools) {
        this.name = name;
        this.id = id;
        this.version = version;
        this.author = author;
        this.description = description;
        this.help = help;
        this.fileParameters = fileParameters;
        this.subTools = subTools;
        this.baseUrl = baseUrl;
    }

    public static class SubToolDescriptor {
        String name;//	Name of the web tool.
        String id;//	 	LSID which uniquely identifies the web tool.
        String version;//	 	The version of the web tool.
        String author;//	 	The author of the web tool.
        String description;//	 	Description of the web tool.
        String help;//	 	URI for help about the web tool.
        List<FileParameter> fileParameters;//		A list of FileParameters.  See below.
        String urlModifier; // 	The URL fragment that will be appended to the parent WebTool’s baseUrl.

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
    }

    public static class FileParameter {
        String name;// 	Name of the parameter.
        String description;// 	Description of the parameter.
        String required;// 	Boolean flag indicating whether or not the parameter is required.
        String compositeFilename;//  	True means that the value for this parameter can contain multiple filenames.  The filenames are separated by the nameDelimeters below.
        String nameDelimiters;// 	The string used to delimit multiple filenames.  This is only used if compositeFilename is true.
        List<GSDataFormat> formats;// 	A set of GSDataFormats, describing the acceptable file formats for this FileParameter.  See below.

        public FileParameter(String name, String description, String required, String compositeFilename,
                             String nameDelimiters, List<GSDataFormat> formats) {
            this.name = name;
            this.description = description;
            this.required = required;
            this.compositeFilename = compositeFilename;
            this.nameDelimiters = nameDelimiters;
            this.formats = formats;
        }
    }

    public static class GSDataFormat {

        String name;//	The name of the format.
        String version;//		The version of the format
        String url;//		An URL to the format definition.

        public GSDataFormat(String name, String version, String url) {
            this.name = name;
            this.version = version;
            this.url = url;
        }

    }
}
