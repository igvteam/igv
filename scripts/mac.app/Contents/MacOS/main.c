//
//  main.c
//  igv-launcher
//
//  Created by Jim Robinson on 9/14/24.
//
//  Launch IGV on MacOS (https://github.com/igvteam/igv).
//
//  Code adapted from https://incenp.org/notes/2023/universal-java-app-on-macos.html
//

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <string.h>
#include <err.h>
#include <mach-o/dyld.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <CoreFoundation/CoreFoundation.h>


/* Remove the last n components of a pathname. The pathname
 * is modified in place. Returns 0 if the requested number
 * of components have been removed, -1 otherwise. */
static int remove_last_component(char *buffer , unsigned n)
{
    char *last_slash = NULL;
    while ( n-- > 0 ) {
        if ( (last_slash = strrchr(buffer, '/')) )
            *last_slash = '\0';
    }
    return last_slash ? 0 : -1;
}

/*
 * Test if a given path exists and if it is a  directory.
 * It returns 1 if given path is directory and  exists
 * otherwise returns 0.
 */
static int isDirectoryExists(const char *path)
{
    struct stat stats;
    stat(path, &stats);

    // Check for file existence
    if (S_ISDIR(stats.st_mode))
        return 1;

    return 0;
}

static int simpleLauncher(int argc, char **argv)
{
    char app_path[PATH_MAX];
    uint32_t path_size = PATH_MAX;
    int ret = 0;

    (void) argc;
    (void) argv;

    /* Get the path to the "Contents" directory. */
    if ( _NSGetExecutablePath(app_path, &path_size) == -1 )
        err(EXIT_FAILURE, "Cannot get application directory");
    if ( remove_last_component(app_path, 2) == -1 )
        errx(EXIT_FAILURE, "Cannot get application directory");

    /* Move to that directory. */
    if ( chdir(app_path) == -1 )
        err(EXIT_FAILURE, "Cannot change current directory");

    char* cmd;
    
    /* See if there is a bundled JDK */
    if(isDirectoryExists("jdk-21")) {
        char jdkPath[2048] = {'\0'};
        snprintf(jdkPath, sizeof(jdkPath), "%s/jdk-21", app_path);
        
        setenv("JAVA_HOME", jdkPath, true);
        
        char javaPath[2048] = {'\0'};
        snprintf(javaPath, sizeof(javaPath), "%s/bin/java", jdkPath);
        cmd = javaPath;

    } else {
    
        cmd = "java";
    }
    
    /* Get user defined extra arguments, if any */
    char extraArgumentsPath[2048] = {'\0'};
    const char *homeDir = getenv("HOME");
    snprintf(extraArgumentsPath, sizeof(extraArgumentsPath), "%s/.igv/java_arguments", homeDir);

    if(access(extraArgumentsPath, F_OK) == 0) {
        char extraArguments[2048] = {'\0'};
        snprintf(extraArguments, sizeof(extraArguments), "@%s", extraArgumentsPath);

        /* Run IGV with user extra arguments */

        
        char* args[] = {
            "java",
            "-showversion",
            "-cp",
            "Java/lib/*",
            "-Xmx8g",
            "@Java/igv.args",
            "-Xdock:name=IGV",
            "-Xdock:icon=Resources/IGV_64.png",
            "-Dsamjdk.snappy.disable=true",
            "-Djava.net.preferIPv4Stack=true",
            "-Dapple.laf.useScreenMenuBar=true",
            "-Djava.net.useSystemProxies=true",
            extraArguments,
            "org.broad.igv.ui.Main",
            NULL};

        ret = execvp(cmd, args);
        return ret;

    } else {
        
        /* Run IGV without user extra arguments */

        char* args[] = {
            "java",
            "-showversion",
            "-cp",
            "Java/lib/*",
            "-Xmx8g",
            "@Java/igv.args",
            "-Xdock:name=IGV",
            "-Xdock:icon=Resources/IGV_64.png",
            "-Dsamjdk.snappy.disable=true",
            "-Djava.net.preferIPv4Stack=true",
            "-Dapple.laf.useScreenMenuBar=true",
            "-Djava.net.useSystemProxies=true",
            "org.broad.igv.ui.Main",
            NULL};
        
        ret = execvp(cmd, args);
        return ret;
    }
}

int main(int argc, char **argv)
{
    return simpleLauncher(argc, argv);
}
