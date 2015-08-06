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
package com.iontorrent.data;

/**
 *
 * @author Chantal Roth
 */
public class ReadInfo {
    String readname;
    int flowposition;
    int flowvalue;
    char base;
    
    public ReadInfo(String readname, int flowposition,int flowvalue, char base) {
        this.readname = readname;
        this.flowposition = flowposition;
        this.flowvalue = flowvalue;
        this.base = base;
    }
    public int getFlowPosition() {
        return flowposition;
    }
    public char getBase() {
        return base;
    }
    public String getReadName() {
        return readname;
    }
    public static String getHeader() {
       return "Read name, flow position, base, flow value";
    }
    public String toCsv() {
        StringBuilder b = new StringBuilder();
        b = b.append(readname).append(",").append(flowposition).append(",").append(base).append(",").append(flowvalue);
        return b.toString();
    }
    public String toString() {
        return toCsv();
    }
}
