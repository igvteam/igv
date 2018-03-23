package org.broad.igv.sam.mutreview;

import com.google.gson.Gson;
import com.google.gson.JsonObject;
import org.broad.igv.batch.CommandExecutor;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ga4gh.GoogleUtils;
import org.broad.igv.ga4gh.OAuthUtils;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.sam.AlignmentCounts;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackClickEvent;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.SnapshotUtilities;

import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.util.Optional;

/**
 * Created by jrobinso on 11/14/17.
 */
public class VariantReviewAction {


    public static void scoreMutationItem(final JComponent component, AlignmentTrack track, final TrackClickEvent te) {


        try {

            JsonObject userProfile = OAuthUtils.getInstance().fetchUserProfile();

            if (userProfile == null) {
                MessageUtils.showMessage("Log in via Google Menu before scoring variant");
                return;
            }

            ReferenceFrame frame = te.getFrame();
            VariantReviewMetadata metadata = new VariantReviewMetadata();
            metadata.userEmail = userProfile.get("email").getAsString();
            metadata.userId = userProfile.get("id").getAsString();
            metadata.chrom = frame.getChrName();
            metadata.pos = (int) te.getChromosomePosition() + 1;
            metadata.windowSize = frame.getCurrentRange().getLength();

            AlignmentCounts alignmentCounts = track.getCoverageTrack().getCounts(frame.getChrName(), (int) te.getChromosomePosition(), frame);
            BaseCounts baseCounts = new BaseCounts(alignmentCounts, (int) te.getChromosomePosition());

            byte[] seq = GenomeManager.getInstance().getCurrentGenome().getSequence(frame.getChrName(), (int) te.getChromosomePosition(), (int) te.getChromosomePosition() + 1);
            byte ref = seq[0];
            metadata.ref = (char) ref;
            int r = -1;
            for (int i = 0; i < baseCounts.bases.length; i++) {
                if (baseCounts.bases[i] == ref) {
                    r = i;
                    break;
                }
            }

            metadata.readDepth = baseCounts.totalCount;
            metadata.refCount = baseCounts.posCounts[r] + baseCounts.negCounts[r];
            metadata.refCountPos = baseCounts.posCounts[r];
            metadata.refCountNeg = baseCounts.negCounts[r];

            metadata.alt = "";
            metadata.altCount = "";
            metadata.altCountNeg = "";
            metadata.altCountPos = "";

            for (int i = 0; i < baseCounts.bases.length; i++) {
                if (baseCounts.bases[i] != ref) {
                    if (metadata.alt.length() > 0) {
                        metadata.alt += ",";
                        metadata.altCount += ",";
                        metadata.altCountPos += ",";
                        metadata.altCountNeg += ",";
                    }
                    metadata.alt += baseCounts.bases[i];
                    metadata.altCount += (baseCounts.negCounts[i] + baseCounts.posCounts[i]);
                    metadata.altCountPos += baseCounts.posCounts[i];
                    metadata.altCountNeg += baseCounts.negCounts[i];
                }
            }

            metadata.delCount = baseCounts.delCount;


            // Create the screenshot
            PreferencesManager.forceDefaults = true;
            
            CommandExecutor cmdExe = new CommandExecutor();
            cmdExe.setSleepInterval("0");

            int chrPosition = (int) (te.getChromosomePosition()) + 1;  // Convert to "1" base coords
            int start = chrPosition - 100;
            int end = chrPosition + 100;

            cmdExe.execute("goto " + te.getFrame().getChrName() + ":" + start + "-" + end);
            cmdExe.execute("sort base " + chrPosition);
            Rectangle clipRect = track.alignmentsRect;
            BufferedImage image = SnapshotUtilities.createBufferedImage(component, clipRect, 1000);
            PreferencesManager.forceDefaults = false;

            Frame parent = IGV.hasInstance() ? IGV.getMainFrame() : null;
            VariantReviewFX.open(parent, image, metadata);

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
