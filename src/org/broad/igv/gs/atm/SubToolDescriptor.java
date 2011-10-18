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
