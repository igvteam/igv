package org.broad.igv.ui.commandbar;

import com.jidesoft.hints.ListDataIntelliHints;
import htsjdk.tribble.NamedFeature;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.ui.action.SearchCommand;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.util.LongRunningTask;

import javax.swing.*;
import javax.swing.text.JTextComponent;
import java.util.List;

/**
 * Created by jrobinso on 7/6/17.
 */
public class SearchTextField extends JTextField {

    static Logger log = LogManager.getLogger(SearchTextField.class);

    public SearchTextField() {

        setToolTipText("Enter a gene or locus, e.f. EGFR,   chr1,   or chr1:100,000-200,000");

        addActionListener(actionevent -> {
            searchByLocus(getText());
        });

        new SearchHints(this);  // This has the side-effect, apparently, of enabling hints
    }


    public void searchByLocus(final String searchText) {
        LongRunningTask.submit((new SearchCommand(FrameManager.getDefaultFrame(), searchText)));
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
                //List<IGVNamedFeature> features = SearchCommand.comprehensiveFeatureSearch(text);
                List<NamedFeature> features = GenomeManager.getInstance().getCurrentGenome().getFeatureDB().getFeaturesList(text, SearchCommand.SEARCH_LIMIT);
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
