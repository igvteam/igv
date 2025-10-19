package org.broad.igv.seg;

import java.util.List;

public record SampleGroup(String label, List<String> samples) {
}
