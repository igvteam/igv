package org.igv.session.autosave;

import org.igv.DirectoryManager;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.session.Session;
import org.igv.session.SessionWriter;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.attribute.BasicFileAttributes;
import java.nio.file.attribute.FileTime;
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
        File autosavedSessionFile = new File(DirectoryManager.getAutosaveDirectory(), "exit_session_autosave.xml");
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
        File autosavedSessionFile = new File(DirectoryManager.getAutosaveDirectory(), "exit_session_autosave.xml");
        // Save the session
        (new SessionWriter()).saveSession(session, autosavedSessionFile);
    }

    /**
     * Gets a File object for each timed autosave session file in the autosave directory and returns an array of those
     * objects
     *
     * @return Array of autosave files
     */
    public static synchronized File[] getTimedSessionAutosaveFiles() {
        File[] autosaveFiles = DirectoryManager.getAutosaveDirectory().listFiles();
        // Remove the exit session autosave file from the list
        return Arrays.stream(autosaveFiles)
                .filter(file -> !file.getName().equals("exit_session_autosave.xml"))
                .toArray(File[]::new);
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
        File[] autosaves = getTimedSessionAutosaveFiles();
        if(autosaves.length > PreferencesManager.getPreferences().getAsInt(Constants.AUTOSAVES_TO_KEEP)) {
            Arrays.sort(autosaves);
            autosaves[0].delete();
        }

        // Save the session
        (new SessionWriter()).saveSession(session, autosavedSessionFile);
    }

    /**
     * Gets a File object for the most recent file in the autosave directory, or an empty optional if the directory is
     * empty
     *
     * @return Optional containing a File object for the session autosave file or empty
     */
    public static synchronized Optional<File> getMostRecentAutosaveFile() throws IOException{
        File[] autosaveFiles = DirectoryManager.getAutosaveDirectory().listFiles();
        if(autosaveFiles.length < 1) {
            return Optional.empty();
        }
        // Check the creation dates of the files so we can figure out which is the most recent
        int latestIndex = 0;
        long latestCreationTime = Files.readAttributes(autosaveFiles[0].toPath(), BasicFileAttributes.class)
                .creationTime().toMillis();
        for(int i = 1; i < autosaveFiles.length; i++) {
            // This file is the new latest if its creationTime is later than the current latest
            long currentCreationTime = Files.readAttributes(autosaveFiles[i].toPath(), BasicFileAttributes.class)
                    .creationTime().toMillis();
            if(currentCreationTime > latestCreationTime) {
                latestCreationTime = currentCreationTime;
                latestIndex = i;
            }
        }
        return Optional.of(autosaveFiles[latestIndex]);
    }
}
