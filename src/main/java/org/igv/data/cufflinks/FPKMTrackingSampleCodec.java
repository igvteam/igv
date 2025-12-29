package org.igv.data.cufflinks;


/**
 * Wrapper class to make it easy to just extract the first sample from
 * an fpkm_tracking file.
 * @author jacob
 * @date 2013-May-28
 */
public class FPKMTrackingSampleCodec extends CufflinksCodec<FPKMSampleValue>   {

    private FPKMTrackingCodec trackingCodec;

    public FPKMTrackingSampleCodec() {
        super(FPKMSampleValue.class, "Plugin");
        this.trackingCodec = new FPKMTrackingCodec(this.path);
    }

    @Override
    protected Object readHeader(String[] tokens) {
        return trackingCodec.readHeader(tokens);
    }

    @Override
    public FPKMSampleValue decode(String line) {
        FPKMValue val = trackingCodec.decode(line);
        if(val != null){
            return val.getSampleValue(0);
        }else{
            return null;
        }
    }

    @Override
    public boolean canDecode(String path) {

        String fn = path.toLowerCase();
        if(fn.endsWith(".gz")) fn = fn.substring(0, fn.length()-3);
        return fn.endsWith("fpkm_tracking");
    }
}
