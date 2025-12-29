package org.broad.igv.event;

/**
 * @author Jim Robinson
 * @date 12/2/11
 */
public record AlignmentTrackEvent(Type type) implements IGVEvent{
    public enum Type {ALLELE_THRESHOLD, RELOAD, REFRESH}
}