package org.igv.ui.util;

/**
 * @author Jim Robinson
 * @date 1/25/12
 */
public class ImageFileTypes {

    /**
     * Snapshot types
     */
    public static enum Type {

        NULL("", ""),
        EPS(".eps", "Encapsulated Postscript Files (*.eps)"),
        PDF(".pdf", "Portable Document FormatFles (*.pdf)"),
        SVG(".svg", "Scalable Vector Graphics Files (*.svg)"),
        PNG(".png", "Portable Network Graphics Files (*.png)"),
        JPEG(".jpeg", "Joint Photographic Experts Group Files (*.jpeg)");
        private String fileExtension;
        private String fileDescription;

        Type(String extension, String description) {
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

    public static Type getImageFileType(String fileExtension) {

        String extension = fileExtension.toLowerCase();
        Type type = Type.NULL;

        if(".jpg".equalsIgnoreCase(fileExtension)){
            type = Type.JPEG;
        }

        for(Type iterType: Type.values()){
            if(type != Type.NULL) break;
            if(iterType.getExtension().equalsIgnoreCase(extension)){
                type = iterType;
            }
        }

        return type;
    }

}
