package org.broad.igv.sam;

/**
 * Created by jrobinso on 2/10/17.
 */
public class ExperimentTypeChangeEvent {

    public final Object source;
    public final AlignmentDataManager.ExperimentType type;

    public ExperimentTypeChangeEvent(Object source, AlignmentDataManager.ExperimentType type) {
        this.source = source;
        this.type = type;
    }
}
