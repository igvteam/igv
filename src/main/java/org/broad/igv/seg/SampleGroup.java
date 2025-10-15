package org.broad.igv.seg;

import java.util.List;

public class SampleGroup {
    private String label;
    private List<String> samples;

    public SampleGroup(String label, List<String> samples) {
        this.label = label;
        this.samples = samples;
    }

    public String getLabel() {
        return label;
    }

    public List<String> getSamples() {
        return samples;
    }
}

