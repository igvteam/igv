package org.igv.event;

import org.igv.feature.genome.Genome;

/**
 * Created by jrobinso on 1/7/17.
 */
public record GenomeChangeEvent(Genome genome) implements IGVEvent {}
