package org.broad.igv.ui.form;

public class FormField {

    String key;
    String label;
    String type;
    String defaultValue;
    String comment;

    public FormField(String key, String label, String type, String defaultValue, String comment) {
        this.key = key;
        this.label = label;
        this.type = type;
        this.defaultValue = defaultValue;
        this.comment = comment;
    }

    public String getKey() {
        return key;
    }

    public String getLabel() {
        return label;
    }

    public String getType() {
        return type;
    }

    public String getDefaultValue() {
        return defaultValue;
    }

    public String getComment() {
        return comment;
    }
}
