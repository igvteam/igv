package org.broad.igv.sam;

/**
 * Created by jrobinso on 2/10/17.
 */
public class ExperimentTypeChangeEvent {
    AlignmentDataManager.ExperimentType type;
    public ExperimentTypeChangeEvent(AlignmentDataManager.ExperimentType type) {
        this.type = type;
    }
}
