package org.igv.sample;

import java.util.List;

public record SampleGroup(String label, List<String> samples) {
}
