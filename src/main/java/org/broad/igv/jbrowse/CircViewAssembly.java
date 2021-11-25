package org.broad.igv.jbrowse;

class CircViewAssembly {
    String id;
    String name;
    CircViewRegion [] chromosomes;

    public CircViewAssembly(String id, String name, CircViewRegion[] regions) {
        this.id = id;
        this.name = name;
        this.chromosomes = regions;
    }

    public String toJson() {
        StringBuffer buf = new StringBuffer();
        buf.append("{");
        buf.append(JsonUtils.toJson("id", id));
        buf.append(",");
        buf.append(JsonUtils.toJson("name", name));
        buf.append(",\"chromosomes\":");
        buf.append("[");
        boolean first = true;
        for(CircViewRegion c : this.chromosomes) {
            if(!first) {
                buf.append(",");
            }
            buf.append(c.toJson());
            first = false;
        }
        buf.append("]");
        buf.append("}");
        return buf.toString();
    }
}

class CircViewRegion {
    String name;
    int bpLength;
    String color;
    public CircViewRegion(String name, int bpLength, String color) {
        this.name = name;
        this.bpLength = bpLength;
        this.color = color;
    }

    public String toJson() {
        StringBuffer buf = new StringBuffer();
        buf.append("{");
        buf.append(JsonUtils.toJson("name", name));
        buf.append(",");
        buf.append(JsonUtils.toJson("color", color));
        buf.append(",");
        buf.append(JsonUtils.toJson("bpLength", bpLength));
        buf.append("}");
        return buf.toString();
    }
}
