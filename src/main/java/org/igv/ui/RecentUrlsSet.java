package org.igv.ui;

import org.igv.ui.action.LoadFilesMenuAction;
import org.igv.ui.action.LoadFromURLMenuAction;
import org.igv.util.ResourceLocator;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Objects;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * An implementation of StackSet which supports {@link ResourceLocator}s which include a path and an index
 * This matches the form used by the {@link LoadFilesMenuAction} and {@link LoadFromURLMenuAction}
 *
 *
 */
public class RecentUrlsSet extends StackSet<ResourceLocator> {
    private static final String INDEX_DELIM = "index:";
    private static final Pattern INDEX_SPLITTER = Pattern.compile("\\s" + INDEX_DELIM);

    public RecentUrlsSet(int maxSize) {
        super(maxSize);
    }

    public RecentUrlsSet(Collection<ResourceLocator> c, int maxSize) {
        super(c, maxSize);
    }

    public String asString(){
        return this.stream()
                .map(RecentUrlsSet::locatorToString)
                .collect(Collectors.joining("|"));
    }

    private static String locatorToString(ResourceLocator locator) {
        StringBuilder builder = new StringBuilder();
        builder.append(locator.getPath());
        if(locator.getIndexPath() != null) {
            builder.append(" " + INDEX_DELIM);
            builder.append(locator.getIndexPath());
        }
        return builder.toString();
    }

    private static ResourceLocator stringToLocator(String locationString) {
        String[] split = INDEX_SPLITTER.split(locationString);

        if (split.length != 1 && split.length != 2){
            return null;
        }
        final ResourceLocator result = new ResourceLocator(split[0].strip());
        if(split.length == 2) {
            result.setIndexPath(split[1].strip());
        }
        return result;
    }

    public static RecentUrlsSet fromString(String urls, int maxLength){
        if(urls == null) {
            return new RecentUrlsSet(maxLength);
        }

        String[] elements = urls.split("\\|");
        List<ResourceLocator> locators = Arrays.stream(elements)
                .filter(Objects::nonNull)
                .filter(elem -> !elem.isBlank())
                .map(RecentUrlsSet::stringToLocator)
                .filter(Objects::nonNull)
                .toList();
        return new RecentUrlsSet(locators, maxLength);
    }
}
