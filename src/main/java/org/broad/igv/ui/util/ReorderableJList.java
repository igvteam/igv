package org.broad.igv.ui.util;

import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.Transferable;
import java.awt.dnd.DnDConstants;
import java.awt.dnd.DragGestureEvent;
import java.awt.dnd.DragGestureListener;
import java.awt.dnd.DragGestureRecognizer;
import java.awt.dnd.DragSource;
import java.awt.dnd.DragSourceDragEvent;
import java.awt.dnd.DragSourceDropEvent;
import java.awt.dnd.DragSourceEvent;
import java.awt.dnd.DragSourceListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.List;

import javax.swing.*;

public class ReorderableJList<T> extends JList {
    DefaultListModel model;

    public ReorderableJList() {

        setTransferHandler(new MyListDropHandler());
        new MyDragListener(this);
    }

    public void setElements(List<T> elements) {
        model = new DefaultListModel();
        setModel(model);
        setDragEnabled(true);
        setDropMode(DropMode.INSERT);

        for (T obj : elements) {
            model.addElement(obj);
        }
    }

    public List<T> getElements() {
        List<T> elementList = new ArrayList<T>();
        Enumeration en = model.elements();
        if (en != null) {
            while (en.hasMoreElements()) {
                elementList.add((T) en.nextElement());
            }
        }
        return elementList;
    }

    class MyDragListener implements DragSourceListener, DragGestureListener {
        ReorderableJList list;

        DragSource ds = new DragSource();

        public MyDragListener(ReorderableJList list) {
            this.list = list;
            DragGestureRecognizer dgr = ds.createDefaultDragGestureRecognizer(list,
                    DnDConstants.ACTION_MOVE, this);

        }

        public void dragGestureRecognized(DragGestureEvent dge) {
            StringSelection transferable = new StringSelection(Integer.toString(list.getSelectedIndex()));
            ds.startDrag(dge, DragSource.DefaultCopyDrop, transferable, this);
        }

        public void dragEnter(DragSourceDragEvent dsde) {
        }

        public void dragExit(DragSourceEvent dse) {
        }

        public void dragOver(DragSourceDragEvent dsde) {
        }

        public void dragDropEnd(DragSourceDropEvent dsde) {
            if (dsde.getDropSuccess()) {
                //System.out.println("Succeeded");
            } else {
                //System.out.println("Failed");
            }
        }

        public void dropActionChanged(DragSourceDragEvent dsde) {
        }
    }

    class MyListDropHandler extends TransferHandler {

        public boolean canImport(TransferHandler.TransferSupport support) {
            if (!support.isDataFlavorSupported(DataFlavor.stringFlavor)) {
                return false;
            }
            JList.DropLocation dl = (JList.DropLocation) support.getDropLocation();
            if (dl.getIndex() == -1) {
                return false;
            } else {
                return true;

            }
        }

        public boolean importData(TransferHandler.TransferSupport support) {
            if (!canImport(support)) {
                return false;
            }

            Transferable transferable = support.getTransferable();
            String indexString;
            try {
                indexString = (String) transferable.getTransferData(DataFlavor.stringFlavor);
            } catch (Exception e) {
                return false;
            }

            int fromIndex = Integer.parseInt(indexString);
            JList.DropLocation dl = (JList.DropLocation) support.getDropLocation();
            int toIndex = dl.getIndex();

            if (fromIndex < toIndex) {
                toIndex--;
            }

            Object obj = model.remove(fromIndex);
            model.add(toIndex, obj);
            invalidate();

            return true;
        }
    }


    // Example


    public static void main(String[] a) {
        List<String> list = Arrays.asList("a", "b", "c");
        JDialog f = new JDialog();
        f.setModal(true);
        final ReorderableJList<String> roList = new ReorderableJList<String>();
        roList.setElements(list);
        f.add(new JScrollPane(roList));
        f.setSize(300, 300);
        f.setVisible(true);

        for (String s : roList.getElements()) {
            System.out.println(s);
        }

    }


}