package org.igv.ui;

import org.junit.Test;

import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

/**
 * Command line argument parsing.  A session named on the command line has to be recognized as a session --
 * loading one as a track fails with "Unknown file type".
 */
public class MainArgsTest {

    @Test
    public void testSessionOnItsOwn() {
        for (String path : new String[]{"test/sessions/json/vcf-session.json", "session.xml"}) {
            Main.IGVArgs args = new Main.IGVArgs(new String[]{path});
            assertEquals(path, args.getSessionFile());
            assertTrue(empty(args.getDataFileStrings()));
        }
    }

    /**
     * A session with a locus after it -- the two argument legacy form.
     */
    @Test
    public void testSessionWithLocus() {
        Main.IGVArgs args = new Main.IGVArgs(new String[]{"session.json", "chr1:100-200"});
        assertEquals("session.json", args.getSessionFile());
        assertEquals("chr1:100-200", args.getLocusString());
    }

    /**
     * A session given with other arguments goes through the multi argument path, which treated everything as a
     * data file -- so even an XML session was loaded as a track there.
     */
    @Test
    public void testSessionAmongSeveralArguments() {
        Main.IGVArgs args = new Main.IGVArgs(new String[]{"-l", "chr1:100-200", "session.json", "reads.bam"});
        assertEquals("session.json", args.getSessionFile());
        assertEquals(List.of("reads.bam"), args.getDataFileStrings());
    }

    /**
     * Data files are still data files, alone and among other arguments.
     */
    @Test
    public void testDataFiles() {
        Main.IGVArgs args = new Main.IGVArgs(new String[]{"reads.bam"});
        assertNull(args.getSessionFile());
        assertEquals(List.of("reads.bam"), args.getDataFileStrings());

        args = new Main.IGVArgs(new String[]{"-l", "chr1:100-200", "a.bam", "b.vcf"});
        assertNull(args.getSessionFile());
        assertEquals(List.of("a.bam", "b.vcf"), args.getDataFileStrings());
    }

    /**
     * Only the first session is taken as the session; anything after it is data, as before.
     */
    @Test
    public void testSecondSessionIsNotTakenAsData() {
        Main.IGVArgs args = new Main.IGVArgs(new String[]{"-l", "chr1:100-200", "a.json", "b.json"});
        assertEquals("a.json", args.getSessionFile());
        assertEquals(List.of("b.json"), args.getDataFileStrings());
    }

    private static boolean empty(List<String> list) {
        return list == null || list.isEmpty();
    }
}
