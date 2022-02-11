package org.broad.igv.ui.commandbar;

import com.jidesoft.hints.ListDataIntelliHints;
import org.broad.igv.logging.*;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.NamedFeature;
import org.broad.igv.ui.GlobalKeyDispatcher;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.action.SearchCommand;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;

import javax.swing.*;
import javax.swing.text.JTextComponent;
import java.awt.event.*;
import java.util.List;

/**
 * Created by jrobinso on 7/6/17.
 */
public class SearchTextField extends JTextField {

    static Logger log = LogManager.getLogger(SearchTextField.class);

    public SearchTextField() {

        setToolTipText("Enter a gene or locus, e.f. EGFR,   chr1,   or chr1:100,000-200,000");
        
        this.addFocusListener(new FocusListener() {
            @Override
            public void focusGained(FocusEvent e) {
                GlobalKeyDispatcher.getInstance().disable();
            }

            @Override
            public void focusLost(FocusEvent e) {
                GlobalKeyDispatcher.getInstance().enable();
            }
        });

        addActionListener(actionevent -> {
            searchByLocus(getText());
            GlobalKeyDispatcher.getInstance().enable();
        });

        new SearchHints(this);  // This has the side-effect, apparently, of enabling hints
    }


    public void searchByLocus(final String searchText) {

        (new SearchCommand(FrameManager.getDefaultFrame(), searchText)).execute();

    }

    private class SearchHints extends ListDataIntelliHints<String> {

        public SearchHints(JTextComponent jTextComponent) {
            super(jTextComponent, new String[]{});
        }

        @Override
        public void acceptHint(Object context) {
            String text = (String) context;
            super.acceptHint(context);
            searchByLocus(text);
        }

        @Override
        public boolean updateHints(Object context) {
            String text = (String) context;
            if (text.length() <= 1) {
                return false;
            } else {
                //TODO Uncomment to use comprehensive feature search, note that it should support partial matches
                //List<NamedFeature> features = SearchCommand.comprehensiveFeatureSearch(text);
                List<NamedFeature> features = FeatureDB.getFeaturesList(text, SearchCommand.SEARCH_LIMIT);
                final List<SearchCommand.SearchResult> results = SearchCommand.getResults(features);
                Object[] list = SearchCommand.getSelectionList(results, false);
                if (list.length >= 1) {
                    this.setListData(list);
                    return true;
                }
            }
            return false;
        }
    }


}
