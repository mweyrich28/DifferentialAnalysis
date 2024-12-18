package org.src;

import augmentedTree.IntervalTree;
import net.sf.samtools.AlignmentBlock;

import java.util.*;

public class Transcript {
    private final String transcriptId;
    private final int start;
    private final int stop;
    private final ArrayList<Exon> exonList = new ArrayList<>();
    private final HashMap<Integer, Exon> exonStarts = new HashMap<>();
    private final HashMap<Integer, Exon> exonEnds = new HashMap<>();
    private final char strand;
    private String transcriptSeq; // patched together using its exons
    public Transcript(String transcriptId, char strand, int start, int stop) {
        this.transcriptId = transcriptId;
        this.strand = strand;
        this.start = start;
        this.stop = stop;
    }

    public void addExon(int start, int end, int pos) {
        Exon exon = new Exon(start, end, pos, end - start + 1);
        exonList.add(exon);
        exonStarts.put(start, exon);
        exonEnds.put(end, exon);
    }

    public String getTranscriptId() {
        return transcriptId;
    }

    public ArrayList<Exon> getExonList() {
        return this.exonList;
    }

    public HashMap<Integer, Exon> getExonEnds() {
        return exonEnds;
    }

    public HashMap<Integer, Exon> getExonStarts() {
        return exonStarts;
    }

    public void reversExonList() {
        Collections.reverse(this.exonList);
    }
}
