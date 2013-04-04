/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.gs.atm;

import com.google.gson.annotations.SerializedName;

import java.util.List;

/**
 * @author Jim Robinson
 * @date Aug 2, 2011
 */
public class WebToolDescriptor {
    private String name;

    @SerializedName("internalId") private String id;

    private String author;
    private String description;
    private String baseUrl;
    private List<FileParameter> fileParameters;

    public WebToolDescriptor(String name, String id, String author, String description, String baseUrl,
                             List<FileParameter> fileParams) {

        this.name = name;
        this.id = id;
        this.author = author;
        this.description = description;
        this.fileParameters = fileParams;
        this.baseUrl = baseUrl;
    }


    public void print() {
        System.out.println();
        System.out.println(name);
        for(FileParameter fp : fileParameters) {
            fp.print();
        }

    }

    public String getName() {
        return name;
    }

    public String getId() {
        return id;
    }

    public String getAuthor() {
        return author;
    }

    public String getDescription() {
        return description;
    }

    public String getBaseUrl() {
        return baseUrl;
    }

    public List<FileParameter> getFileParameters() {
        return fileParameters;
    }
}
