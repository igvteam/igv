package org.broad.igv.ui;

/**
 * Backward compatibility shim.  IGV's main class moved to {@link org.igv.ui.Main}; this class remains for
 * launchers, scripts, and manifests that still reference the old package.
 */
public class Main {

    public static void main(final String args[]) {
        org.igv.ui.Main.main(args);
    }

}
