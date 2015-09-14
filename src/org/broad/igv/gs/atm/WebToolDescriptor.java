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
