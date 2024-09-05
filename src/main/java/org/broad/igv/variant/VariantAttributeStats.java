package org.broad.igv.variant;

import org.broad.igv.track.AttributeManager;

import java.util.*;

public class VariantAttributeStats {

    private VariantAttributeStats(){}

    private final Map<Key, Stats> attributes = new TreeMap<>();

    public enum Type{
        INFO, FORMAT
    }

    private record Key(Type type, String name) implements Comparable<Key> {


        final static Comparator<Key> COMPARATOR = Comparator.comparing(Key::type)
                .thenComparing(Key::name);

        @Override
        public int compareTo(Key o) {
            return COMPARATOR.compare(this, o);
        }
    };

    public sealed interface Stats {
        void add(Object value);

        final class Discrete implements Stats{
            private final Map<String, Integer> values = new HashMap<>();
            private static final int MAX_SIZE = 20;

            @Override
            public void add(Object value) {
                if( MAX_SIZE < values.size()){
                    values.compute(value.toString(), (k, v) -> v == null ? 1 : v + 1);
                }
            }

            public List<String> getValues(){
                return values.entrySet().stream()
                        .sorted(Comparator.comparing(Map.Entry::getValue))
                        .map(Map.Entry::getKey)
                        .toList();
            }

        }
        final class Continuous implements Stats{
            @Override
            public void add(Object value) {
                //
            }
        };
    }

    public void add(Variant v){
       for( Map.Entry<String, Object> entry: v.getAttributes().entrySet()){
           //attributes.compute()//
        }
    }



    public Stats getStats(Type type, String name) {return attributes.get(new Key(type, name));}
    private static VariantAttributeStats instance;

    public static VariantAttributeStats getInstance() {
        if (instance == null) {
            instance = new VariantAttributeStats();
        }
        return instance;
    }


}
