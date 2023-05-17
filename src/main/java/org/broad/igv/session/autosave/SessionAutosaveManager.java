package org.broad.igv.session.autosave;

import org.broad.igv.DirectoryManager;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.session.Session;
import org.broad.igv.session.SessionWriter;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Defines several static methods for saving and retrieving session autosave files
 */
public class SessionAutosaveManager {


    /**
     * Returns a file object for the file corresponding to the autosave-on-exit session file, or an empty optional if
     * the file does not exist
     *
     * @return Optional containing a File object for the session autosave file or empty
     */
    public static synchronized Optional<File> getExitSessionAutosaveFile() {
        File autosavedSessionFile = new File(DirectoryManager.getAutosaveDirectory(), "session_autosave.xml");
        if(!autosavedSessionFile.exists()) {
            return Optional.empty();
        }
        return Optional.of(autosavedSessionFile);
    }

    /**
     * Saves the specified session to a file named to indicate it is an autosave.  Creates the file if one does not yet
     * exist, or overwrites an existing session autosave created on exit
     *
     * @param session the session to save
     * @throws IOException if saving the session fails
     */
    public static synchronized void saveExitSessionAutosaveFile(Session session) throws IOException {
        // Get the file we use for saving sessions on exit
        File autosavedSessionFile = new File(DirectoryManager.getAutosaveDirectory(), "session_autosave.xml");
        // Save the session
        (new SessionWriter()).saveSession(session, autosavedSessionFile);
    }

    /**
     * Gets a File object for each file in the autosave directory and returns an array of those objects
     *
     * @return Array of autosave files
     */
    public static synchronized File[] getSessionAutosaveFiles() {
        return DirectoryManager.getAutosaveDirectory().listFiles();
    }

    /**
     * Creates a new timed autosave session file and writes the provided session to it
     *
     * @param session The session to save
     * @throws IOException if saving the session or deleting sessions over the limit fails
     */
    public static synchronized void saveTimedSessionAutosaveFile(Session session) throws IOException {
        // Create a filename for the new autosave file that includes the current datetime
        String currentDate = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss.SSS'Z'").format(new Date());
        String newAutosaveFilename = "session_autosave" + currentDate + ".xml";

        // Create the file
        File autosavedSessionFile = new File(DirectoryManager.getAutosaveDirectory(), newAutosaveFilename);
        autosavedSessionFile.createNewFile();

        // If this puts the number of timed autosave files over the maximum, delete the oldest (the first in
        // alphabetical order that is the one created on exit)
        List<File> autosaves = Arrays.stream(DirectoryManager.getAutosaveDirectory().listFiles())
                .filter(file -> !file.getName().equals("session_autosave.xml"))
                .collect(Collectors.toList());
        if(autosaves.size() > PreferencesManager.getPreferences().getAsInt(Constants.AUTOSAVES_TO_KEEP)) {
            Collections.sort(autosaves);
            autosaves.get(0).delete();
        }

        // Save the session
        (new SessionWriter()).saveSession(session, autosavedSessionFile);
    }
}
