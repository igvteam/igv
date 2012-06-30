/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package com.iontorrent.data;

import java.util.ArrayList;
import java.util.TreeMap;

/**
 *
 * @author Chantal Roth
 */
public class FlowDistribution {

    private TreeMap<Short, Integer> map;
    private String information;
    private String name;
    private char base;
    private int nrflows;
    /**
     * the chromosome location
     */
    private int location;
    private ArrayList<ReadInfo> readinfos;
    private boolean forward;
    private boolean reverse;

    public FlowDistribution(int location, int nrflows, TreeMap<Short, Integer> map, String name, char base, boolean forward, boolean reverse, String information) {
        this.map = map;
        this.information = information;
        this.name = name;
        this.nrflows = nrflows;
        this.location = location;
        this.forward = forward;
        this.reverse = reverse;
        this.base = base;
    }

    public int getNrFlows() {
        return nrflows;
    }

    public String getName() {
        return name;
    }

    public String toCsv(int binsize) {
        int[] bins = getBinnedData(binsize);
        String nl = "\n";
        StringBuilder csv = new StringBuilder();
        csv = csv.append(getInformation());
        csv = csv.append(nl).append("flow value, count\n");
        for (int b = 0; b < bins.length; b++) {
            csv = csv.append(b * binsize).append(",").append(bins[b]).append(nl);
        }
        csv = csv.append(nl);
        return csv.toString();
    }

    public String toJson() {
        StringBuilder buf = new StringBuilder();
        buf.append("{\n");
        for (Short key : map.keySet()) {
            buf.append("    \"").append(key).append("\" : \"").append(map.get(key)).append("\"\n");
        }
        buf.append("}\n");
        return buf.toString();
    }

    public String getReadInfoString() {
        String nl = "\n";
        StringBuilder csv = new StringBuilder();
        csv = csv.append(getInformation());
        csv = csv.append(nl).append(ReadInfo.getHeader()).append(nl);
        for (ReadInfo ri : readinfos) {
            csv = csv.append(ri.toCsv()).append(nl);
        }
        csv = csv.append(nl);
        return csv.toString();
    }
    public String getReadNames() {
        StringBuilder names = new StringBuilder();
        for (ReadInfo ri : readinfos) {
            names = names.append(ri.getReadName()).append("_");
        }
      
        return names.toString();
    }

    public int[] getBinnedData(int binsize) {
        int maxx = 0;
        for (Short x : map.keySet()) {
            if (x > maxx) {
                maxx = x;
            }
        }
        int nrbins = maxx / binsize + 1;
        int bins[] = new int[nrbins];
        for (Short x : map.keySet()) {
            int y = map.get(x);
            bins[x / binsize] += y;
        }
        return bins;
    }

    /**
     * @return the map
     */
    public TreeMap<Short, Integer> getMap() {
        return map;
    }

    /**
     * @param map the map to set
     */
    public void setMap(TreeMap<Short, Integer> map) {
        this.map = map;
    }

    /**
     * @return the name
     */
    public String getInformation() {
        return information;
    }

    /**
     * @return the location
     */
    public int getLocation() {
        return location;
    }

    /**
     * @param location the location to set
     */
    public void setLocation(int location) {
        this.location = location;
    }

    public int getMaxX() {
        int maxx = 0;
        for (Short x : map.keySet()) {
            if (x > maxx) {
                maxx = x;
            }
        }
        return maxx;
    }

    public void setReadInfos(ArrayList<ReadInfo> readinfos) {
        this.readinfos = readinfos;
    }

    public char getBase() {
        return base;
    }
    public boolean isForward() {
        return forward;
    }
    public boolean isReverse() {
        return reverse;
    }
}
