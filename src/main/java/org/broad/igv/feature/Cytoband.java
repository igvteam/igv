package org.broad.igv.feature;


public class Cytoband implements IGVNamedFeature {
    String chromosome;
    String name;
    int end;
    int start;
    char type; // p, n, or c
    short stain;


    public Cytoband(String chromosome) {
        this.chromosome = chromosome;
        this.name = "";
    }

    public Cytoband(String chromosome, int start,int end, String name, String gieStain) {
        this.chromosome = chromosome;
        this.end = end;
        this.start = start;
        this.name = name;
        if (gieStain.equals("acen")) {
            setType('c');
        } else {
            setType(gieStain.charAt(1));
            if (type == 'p') {
                String stainString = gieStain.substring(4).trim();
                short stain = stainString.length() == 0 ? 100 : Short.parseShort(stainString);
                setStain(stain);
            }
        }
    }

    public void trim() {

        // @todo -- trim arrays
    }

    @Override
    public String getContig() {
        return chromosome;
    }

    public String getChr() {
        return chromosome;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    public int getEnd() {
        return end;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public int getStart() {
        return start;
    }

    public void setType(char type) {
        this.type = type;
    }

    public char getType() {
        return type;
    }

    public void setStain(short stain) {
        this.stain = stain;
    }

    public short getStain() {
        return stain;
    }


}

