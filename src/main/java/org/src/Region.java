package org.src;

import augmentedTree.Interval;

public class Region implements Interval {
    private int start;
    private int stop;
    private String id;
    private String transcriptId;
    private String type = null;
    private int pos;

    public Region(int start, int end) {
        this.start = start;
        this.stop = end;
    }
    public Region(int start, int end, int pos) {
        this.start = start;
        this.stop = end;
        this.pos = pos;
    }

    public Region(String id, int start, int end) {
        this.id = id;
        this.start = start;
        this.stop = end;
    }

    @Override
    public int getStart() {
        return start;
    }

    @Override
    public int getStop() {
        return stop;
    }

    public void setStop(int stop) {
        this.stop = stop;
    }

    public void setStart(int start) {
        this.start = start;
    }

    @Override
    public String toString() {
        return getId() + ": " + getStart() + "-" + getStop() + " → " + getTranscriptId() + " → " + getType();
    }
    public int getPos() {
        return pos;
    }

    public void setPos(int pos) {
        this.pos = pos;
    }

    @Override
    public int hashCode() {
        int result = 17;
        result = 31 * result + start;
        result = 31 * result + stop;
        result = 31 * result + (id != null ? id.hashCode() : 0);
        return result;
    }
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        Region region = (Region) obj;
        return start == region.start && stop == region.stop;
    }
    public boolean contains(Region other) {
        return this.start <= other.start && this.stop >= other.stop;
    }

    public String getId() {
        return id;
    }

    public void setTranscriptId(String transcriptId) {
        this.transcriptId = transcriptId;
    }

    public String getTranscriptId() {
        return transcriptId;
    }

    public String getType() {
        return this.type;
    }
    public void setType(String t) {
        this.type = t;
    }

    public String debug() {
        return getId() + ": " + getStart() + "-" + getStop() + " → " + getTranscriptId();
    }
}
