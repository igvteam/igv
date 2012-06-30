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
