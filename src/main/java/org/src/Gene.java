package org.src;

import augmentedTree.Interval;
import augmentedTree.IntervalTree;
import net.sf.samtools.SAMRecord;

import java.util.*;

public class Gene implements Interval {
    private final int start;
    private final int end;
    private int meltedLength = 0;

    private final String geneId;
    private final ArrayList<Transcript> transcriptList;
    private final String geneName;
    private final String bioType;
    private final String chr;
    private final char strand;
    private TreeSet<Region> meltedRegions;
    private HashSet<SAMRecord> mappedReads = new HashSet<>();
    private final ArrayList<Intron> introns = new ArrayList<>();
    private IntervalTree<Region> mappedAliBlocksTree = null;
//    private TreeSet<Region> mappedAliBlockSet = new TreeSet<>(Comparator.comparingInt(Region::getStart).thenComparingInt(Region::getStop));
    private ArrayList<Region> mappedAliBlockSet = new ArrayList<>();
    private final IntervalTree<Region> gappedAliBlocksTree = new IntervalTree<>();
    public Gene(String geneId, int start, int end, String geneName, String chr, char strand, String bioType) {
        this.geneId = geneId;
        this.geneName = geneName;
        this.chr = chr;
        this.start = start;
        this.end = end;
        this.strand = strand;
        this.bioType = bioType;
        this.transcriptList = new ArrayList<>();
    }

    public String getGeneId() {
        return geneId;
    }

    public void generateIntrons() {
        for (Transcript transcript : transcriptList) {
            for (int i = 0; i < transcript.getExonList().size() - 1; i++) {
                int intronStart = transcript.getExonList().get(i).getStop() + 1;
                int intronEnd = transcript.getExonList().get(i + 1).getStart() - 1;
                Intron intron = new Intron(intronStart, intronEnd);
                introns.add(intron);
            }
        }
    }

    public void invertTranscripts() {
        for (int i = 0; i < transcriptList.size(); i++) {
            Transcript currTranscript = transcriptList.get(i);
            currTranscript.reversExonList();
            for (int j = 0; j < currTranscript.getExonList().size(); j++) {
                currTranscript.getExonList().get(j).setPos(j);
            }
        }
    }

    public void addTranscript(Transcript transcript) {
        transcriptList.add(transcript);
    }

    public ArrayList<Transcript> getTranscriptList() {
        return transcriptList;
    }

    public String getChr() {
        return chr;
    }

    @Override
    public int getStart() {
        return start;
    }

    @Override
    public int getStop() {
        return end;
    }

    public int getEnd() {
        return end;
    }
    public Transcript getLastTranscript() {
        if (!transcriptList.isEmpty()) {
            return transcriptList.get(transcriptList.size() - 1);
        }
        return null;
    }

    public char getStrand() {
        return strand;
    }

    public String getBioType() {
        return bioType;
    }

    public ArrayList<Exon> getSkippedExons() {
        // store events here
        ArrayList<Exon> skippedExons = new ArrayList<>();
        boolean atLeastOneWT = false;
        // for every intron check every transcript
        for (Intron intron : introns) {
            int intronStart = intron.getStart();
            int intronEnd = intron.getEnd();

            for (Transcript currTranscript : transcriptList) {
                // get relevant HashMaps and check if currTranscript has cds starting or ending at i_s i_e
                HashMap<Integer, Exon> exonEnds = currTranscript.getExonEnds();
                HashMap<Integer, Exon> exonStarts = currTranscript.getExonStarts();

                boolean hasExonInFront = exonEnds.containsKey(intronStart - 1);
                boolean hasExonBehind = exonStarts.containsKey(intronEnd + 1);

                if (hasExonInFront && hasExonBehind) {
                    Exon frontExon = exonEnds.get(intronStart - 1);
                    Exon behindExon = exonStarts.get(intronEnd + 1);

                    // get offset / look if there are exons inbetween
                    int offset = behindExon.getPos() - frontExon.getPos();

                    if (offset != 1) {
                        atLeastOneWT = true;
                        ArrayList<Exon> exonList = currTranscript.getExonList();
                        for (int i = frontExon.getPos(); i < behindExon.getPos(); i++) {
                            if (i > frontExon.getPos() && i < behindExon.getPos()) {
                                skippedExons.add(exonList.get(i));
                            }
                        }
                    }
                }
            }
        }
        if (atLeastOneWT) {
            return skippedExons;
        }
        return null;
    }


    public void melt() {

        ArrayList<Exon> allExons = new ArrayList<>();
        for (Transcript transcript: transcriptList) {
            allExons.addAll(transcript.getExonList());
        }

        Collections.sort(allExons, Comparator.comparingInt(Exon::getStart));
        TreeSet<Region> meltedRegions = new TreeSet<>(
                Comparator.comparingInt(Region::getStart)
                        .thenComparingInt(Region::getStop)
        );


        if (!allExons.isEmpty()) {
            Exon first = allExons.getFirst();
            Region current = new Region(first.getStart(), first.getStop());

            for (int i = 1; i < allExons.size(); i++) {
                Exon exon = allExons.get(i);

                if (exon.getStart() <= current.getStop() + 1) {
                    current.setStop(Math.max(current.getStop(), exon.getStop()));
                } else {
                    meltedRegions.add(current);
                    current = new Region(exon.getStart(), exon.getStop());
                }
            }

            meltedRegions.add(current);
        }
        this.meltedRegions = meltedRegions;
    }

    public TreeSet<Region> getMeltedRegions() {
        if (this.meltedRegions == null) {
            melt();
        }
        return this.meltedRegions;
    }

    public int getMeltedLength() {
        // make sure to only calculate once
        if (this.meltedLength == 0) {
            TreeSet<Region> meltedRegions = getMeltedRegions();
            int length = 0;
            for (Region region : meltedRegions) {
                length += region.getStop() - region.getStart() - 1;
            }
            this.meltedLength = length;
        }
        return meltedLength;
    }

    public void addRead (SAMRecord read) {
        this.mappedReads.add(read);

    }

    public void addAlignedBlocks(TreeSet<Region> blocks) {
        if (blocks.size() > 1) {
            int count = 1;
            int gapStart = 0;
            int gapEnd;
            for (Region r : blocks) {
                if (count % 2 != 0) {
                    gapStart = r.getStop() + 1;
                    count++;
                } else {
                    gapEnd = r.getStart() - 1;
                    Region gap = new Region(r.getId(), gapStart, gapEnd);
                    gappedAliBlocksTree.add(gap);
                    count = 1;
                }

                this.mappedAliBlockSet.add(new Region(r.getId(), r.getStart(), r.getStop()));
            }
            // add last gap if block had uneven length
            if (count == 2) {
                Region r = blocks.getLast();
                gapEnd = r.getStart() - 1;
                Region gap = new Region(r.getId(), gapStart, gapEnd);
                gappedAliBlocksTree.add(gap);
            }
        }
        else { // theres only one alignment block and no gaps
            Region a = blocks.getLast();
            this.mappedAliBlockSet.add(new Region(a.getId(), a.getStart(), a.getStop()));
        }
    }

    public IntervalTree<Region> getMappedAliBlocksTree() {
        if (this.mappedAliBlocksTree == null) {
            this.mappedAliBlocksTree = new IntervalTree<>();
            for (Region r : mappedAliBlockSet) {
                mappedAliBlocksTree.add(r);
            }
            this.mappedAliBlockSet.clear();
        }
        return this.mappedAliBlocksTree;
    }

    public IntervalTree<Region> getGappedAliBlocksTree() {
        return gappedAliBlocksTree;
    }

    public HashSet<SAMRecord> getMappedReads() {
        return mappedReads;
    }
}

