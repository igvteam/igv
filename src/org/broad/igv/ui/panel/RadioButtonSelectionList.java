/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
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
