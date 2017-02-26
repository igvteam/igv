package org.broad.igv.sam;

/**
 * Created by jrobinso on 2/10/17.
 */
public class ExperimentTypeChangeEvent {

    public final AlignmentDataManager.ExperimentType oldType;
    public final AlignmentDataManager.ExperimentType type;

    public ExperimentTypeChangeEvent(AlignmentDataManager.ExperimentType oldType, AlignmentDataManager.ExperimentType type) {
        this.oldType = type;
        this.type = type;
    }
}
