package org.src;

import augmentedTree.Interval;

public class Exon implements Interval {
    private int length;
    private final int start;
    private final int stop;
    private int pos;
    private int IRC;
    private int ERC;
    private String transcriptId;

    public Exon(int start, int end, int pos, int length) {
        this.length = length;
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
        return this.start + "-" + this.stop + " " + "[" + this.pos +"] Length:" + this.length;
    }
    public int getLength() {
        return length;
    }
}
