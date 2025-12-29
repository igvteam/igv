package org.igv.prefs;

import javax.swing.*;

public class PreferencesTextField<T extends JTextField>{
    T obj;
    PreferencesTextField(T obj) {
        this.obj = obj;
    }
    public String getPreferenceText() {
        if (obj.getClass() == JTextField.class) {
            JTextField f = (JTextField) obj;
            return  f.getText();
        }
        if (obj.getClass() == JPasswordField.class) {
            JPasswordField f = (JPasswordField) obj;
            return String.valueOf(f.getPassword());
        }
        return "";
    }

    public T get() {
        return this.obj;
    }
}
