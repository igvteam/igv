/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.ui.panel;

import com.jidesoft.swing.CheckBoxList;
import com.jidesoft.swing.CheckBoxListCellRenderer;
import com.jidesoft.swing.NullRadioButton;

import javax.swing.*;
import java.awt.*;

/** Just like a CheckBoxList, but renders with radio buttons
 * and only allows single selection.
 * User: jacob
 * Date: 2012-Jul-17
 */
public class RadioButtonSelectionList extends CheckBoxList {

    public RadioButtonSelectionList() {
        super();
    }

    @Override
    protected void init() {
        super.init();
        getCheckBoxListSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    }

    @Override
    protected CheckBoxListCellRenderer createCellRenderer() {
        return new RadioButtonListCellRenderer();
    }

    @Override
    protected Handler createHandler() {
        return new SingleSelectionHandler(this);
    }


    protected static class SingleSelectionHandler extends CheckBoxList.Handler {

        public SingleSelectionHandler(CheckBoxList list) {
            super(list);
        }

        protected void toggleSelection(int index) {
            if (_list.getCheckBoxListSelectedIndex() == index) {
                return;
            }

            super.toggleSelection(index);

            //Once we check a box, want to select it also
            _list.setSelectedIndex(_list.getCheckBoxListSelectedIndex());
        }
    }

    public static class RadioButtonListCellRenderer extends CheckBoxListCellRenderer {

        protected AbstractButton button = new NullRadioButton();

        public RadioButtonListCellRenderer() {
            this(null);
        }

        public RadioButtonListCellRenderer(ListCellRenderer renderer) {
            super(renderer);
            //Really have no idea why this is necessary, or how the check box gets added in the first place
            if (getComponentCount() > 0)
                remove(0);
            button.setBorder(BorderFactory.createEmptyBorder(0, 2, 0, 2));
            button.setOpaque(false);
            add(button, BorderLayout.BEFORE_LINE_BEGINS);
            set_checkBox(button);
        }

        public void set_checkBox(AbstractButton button){
            //This is here for JFormDesigner
            //It gets a NoSuchFieldError when trying to access _checkBox
            //Doesn't really matter to the compiled program, only makes
            //editing the dialog easier
            try{
                this._checkBox = button;
            }catch(NoSuchFieldError e){

            }

        }
    }
}
