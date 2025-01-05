package org.src;

import augmentedTree.Interval;

public class Exon implements Interval {
    private int length;
    private final int start;
    private final int stop;
    private int pos;
    private String transcriptId;
    private String oriID; // for debugging

    public Exon(int start, int end, int pos, int length, String ori) {
        this.length = length;
        this.oriID = ori;
        this.start = start;
        this.stop = end;
        this.pos = pos;
    }

    public int getStart() {
        return start;
    }

    public int getStop() {
        return stop;
    }

    public int getPos() {
        return pos;
    }

    public void setPos(int pos) {
        this.pos = pos;
    }

    @Override
    public String toString() {
        return this.oriID + ": " + this.start + "-" + this.stop + " " + "[" + this.pos +"] Length:" + this.length + " → " + getTranscriptId();
    }
    public int getLength() {
        return length;
    }

    public String getTranscriptId() {
        return transcriptId;
    }

    public void setTranscriptId(String transcriptId) {
        this.transcriptId = transcriptId;
    }
}
