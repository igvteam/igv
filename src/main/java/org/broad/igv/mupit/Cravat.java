package org.broad.igv.mupit;

import com.google.gson.Gson;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import javafx.application.Platform;
import javafx.embed.swing.JFXPanel;
import javafx.geometry.Pos;
import javafx.scene.Node;
import javafx.scene.Scene;
import javafx.scene.control.*;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.scene.text.Text;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferenceEditorFX;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.util.BrowserLauncher;
import org.broad.igv.util.HttpUtils;

import javax.swing.*;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by jrobinso on 2/17/17.
 * <p>
 * Static methods for interacting with the "Cravat" webservice
 * <p>
 * http://www.cravat.us/CRAVAT/rest/service/query?mutation=chr22_30421786_+_A_T
 */
public class Cravat {

    public static void main(String[] args) throws IOException {
        test();
    }

    static Color lightGray = new Color(0.9, 0.9, 0.9, 0.5);

    static void test() throws IOException {

        String jsonString = HttpUtils.getInstance().getContentsAsJSON(
                HttpUtils.createURL("http://www.cravat.us/CRAVAT/rest/service/query?mutation=chr22_30421786_+_A_T"));

        JsonParser parser = new JsonParser();
        JsonObject obj = parser.parse(jsonString).getAsJsonObject();

        openCravatView(obj);

    }


    static void openCravatView(JsonObject jsonObject) {


        SwingUtilities.invokeLater(() -> {
            JFrame frame = new JFrame("Cravat");
            final JFXPanel fxPanel = new JFXPanel();
            frame.add(fxPanel);
            frame.setSize(800, 800);
            frame.setVisible(true);
            frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

            Platform.runLater(() -> initFX(fxPanel, jsonObject));
        });

    }


    private static void initFX(JFXPanel fxPanel, JsonObject jsonObject) {


        GridPane gridPane = new GridPane();
        gridPane.setHgap(5);
        gridPane.setVgap(5);
        ScrollPane scrollPane = new ScrollPane(gridPane);
        Scene scene = new Scene(scrollPane);
        int row = 1;
        for (Map.Entry<String, JsonElement> entry : jsonObject.entrySet()) {

            String key = entry.getKey();
            String value = entry.getValue().getAsString();

            final Label keyLabel = new Label(key);
            Node valueLabel;
            if("dbSNP".equals(key)) {
                String link = "https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?searchType=adhoc_search&type=rs&rs=" + value;
                valueLabel = new Hyperlink(value);
                ((Hyperlink) valueLabel).setOnAction(event -> {
                    try {
                        BrowserLauncher.openURL(link);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                });
            }
            else{
                valueLabel = new Label(value);
            }

            StackPane keyPane = new StackPane(keyLabel);
            keyPane.setAlignment(Pos.CENTER_LEFT);
            StackPane valuePane = new StackPane(valueLabel);
            valuePane.setAlignment(Pos.CENTER_LEFT);

            if(row % 2 == 0) {
                keyPane.setStyle("-fx-background-color: #FFFFFF;");
                valuePane.setStyle("-fx-background-color: #FFFFFF;");
            }

            gridPane.add(keyPane, 1, row);
            gridPane.add(valuePane, 2, row);

            row++;

        }

        fxPanel.setScene(scene);
    }
}
