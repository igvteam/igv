package org.broad.igv.synteny;

public class Anchor extends AbstractMapping {

    boolean psuedo;

    public double mapPosition(int fromPosition) {
        double delta = this.scaleFactor * (fromPosition - this.fromStart);

        if (this.direction == true) {
            return this.toStart + delta;
        }

        return this.toEnd - delta;
    }
}