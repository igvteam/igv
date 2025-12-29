package org.broad.igv.util;

/**
 * Enum for indicating how a long running computation completed.
 *
 * @author jacob
 * @date 2013-Oct-23
 */
public enum RunnableResult {
    SUCCESS,
    FAILURE,
    CANCELLED;

    public boolean isSuccess() {
        return this == SUCCESS;
    }
}
