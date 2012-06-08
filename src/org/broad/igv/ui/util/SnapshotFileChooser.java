package org.broad.igv.ui.util;

import org.broad.igv.util.Utilities;

import javax.swing.*;
import javax.swing.filechooser.FileFilter;
import java.awt.*;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.util.LinkedHashMap;

/**
 * @author Jim Robinson
 * @date 1/25/12
 */
public class SnapshotFileChooser extends JFileChooser {


    private static LinkedHashMap<SnapshotFileType, SnapshotFileFilter> SNAPSHOT_TYPE_TO_FILTER = new LinkedHashMap();

    static {

        SNAPSHOT_TYPE_TO_FILTER.put(SnapshotFileType.JPEG, new SnapshotFileFilter(SnapshotFileType.JPEG));
        //SNAPSHOT_TYPE_TO_FILTER.put(SnapshotFileType.PDF,
        //    new SnapshotFileFilter(SnapshotFileType.PDF));
        //SNAPSHOT_TYPE_TO_FILTER.put(SnapshotFileType.EPS,
        //        new SnapshotFileFilter(SnapshotFileType.EPS));
        SNAPSHOT_TYPE_TO_FILTER.put(SnapshotFileType.SVG, new SnapshotFileFilter(SnapshotFileType.SVG));
        SNAPSHOT_TYPE_TO_FILTER.put(SnapshotFileType.PNG, new SnapshotFileFilter(SnapshotFileType.PNG));
    }


    /**
     * Snapshot types
     */
    public static enum SnapshotFileType {

        NULL("", ""),
        EPS(".eps", "Encapsulated Postscript Files (*.eps)"),
        PDF(".pdf", "Portable Document FormatFles (*.pdf)"),
        SVG(".svg", "Scalable Vector Graphics Files (*.svg)"),
        PNG(".png", "Portable Network Graphics Files (*.png)"),
        JPEG(".jpeg", "Joint Photographic Experts Group Files (*.jpeg)");
        private String fileExtension;
        private String fileDescription;

        SnapshotFileType(String extension, String description) {
            fileExtension = extension;
            fileDescription = description;
        }

        public String getExtension() {
            return fileExtension;
        }

        public String getDescription() {
            return fileDescription;
        }
    }

    public static SnapshotFileType getSnapshotFileType(String fileExtension) {

        String extension = fileExtension.toLowerCase();
        SnapshotFileType type = SnapshotFileType.NULL;

        if(".jpg".equalsIgnoreCase(fileExtension)){
            type = SnapshotFileType.JPEG;
        }

        for(SnapshotFileType iterType: SnapshotFileType.values()){
            if(type != SnapshotFileType.NULL) break;
            if(iterType.getExtension().equalsIgnoreCase(extension)){
                type = iterType;
            }
        }

        return type;
    }


    /**
     * Snapshot file filter
     */
    public static class SnapshotFileFilter extends FileFilter {

        private SnapshotFileType type = SnapshotFileType.EPS;

        public SnapshotFileFilter(SnapshotFileType type) {
            this.type = type;
        }

        public boolean accept(File file) {
            if (file.isDirectory()) {
                return true;
            }
            return file.getName().toLowerCase().endsWith(type.getExtension());
        }

        public String getDescription() {
            return type.getDescription();
        }

        public String getExtension() {
            return type.getExtension();
        }

        public boolean accept(File file, String name) {
            return name.toLowerCase().endsWith(type.getExtension());
        }
    }



    boolean accepted = false;
    File previousFile;

    public SnapshotFileChooser(File directory, File selectedFile) {
        super(directory);
        setPreviousFile(selectedFile);
        init();
    }


    public void approveSelection() {
        accepted = true;
        super.approveSelection();
    }

    public void setPreviousFile(File file) {
        this.previousFile = file;
        setSelectedFile(previousFile);
    }

    public File getPreviousFile() {
        return previousFile;
    }

    @Override
    public void cancelSelection() {
        setSelectedFile(null);
        super.cancelSelection();
    }


    @Override
    protected JDialog createDialog(Component parent) throws HeadlessException {
        JDialog dialog = super.createDialog(parent);
        dialog.setLocation(300, 200);
        dialog.setResizable(false);
        dialog.addWindowListener(new WindowAdapter() {

            @Override
            public void windowClosing(WindowEvent e) {
                if (!accepted) {
                    setSelectedFile(null);
                }
            }
        });
        return dialog;
    }

    private void init() {

        FileFilter[] fileFilters =
                SNAPSHOT_TYPE_TO_FILTER.values().toArray(new FileFilter[SNAPSHOT_TYPE_TO_FILTER.size()]);
        // Setup FileFilters
        if (fileFilters != null) {
            for (FileFilter fileFilter : fileFilters) {
                addChoosableFileFilter(fileFilter);
            }
        }

        addPropertyChangeListener(new PropertyChangeListener() {
            public void propertyChange(PropertyChangeEvent e) {


                File oldFile = null;
                String property = e.getPropertyName();
                if (JFileChooser.SELECTED_FILE_CHANGED_PROPERTY.equals(property)) {
                    oldFile = (File) e.getOldValue();
                } else if (JFileChooser.FILE_FILTER_CHANGED_PROPERTY.equals(property)) {

                    if (e.getOldValue() instanceof SnapshotFileFilter &&
                            e.getNewValue() instanceof SnapshotFileFilter) {

                        SnapshotFileFilter newFilter = (SnapshotFileFilter) e.getNewValue();

                        File currentDirectory = getCurrentDirectory();
                        File previousFile = getPreviousFile();
                        if (previousFile != null) {

                            File file = null;
                            if (currentDirectory != null) {
                                file = new File(currentDirectory, previousFile.getName());
                            } else {
                                file = previousFile;
                            }

                            final File selectedFile = Utilities.changeFileExtension(
                                    file, newFilter.getExtension());

                            UIUtilities.invokeOnEventThread(new Runnable() {

                                public void run() {
                                    setPreviousFile(selectedFile);
                                    validate();
                                }
                            });
                        }

                    }
                }
            }
        });
    }
}
