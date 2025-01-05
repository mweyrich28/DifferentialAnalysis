package org.src;

import augmentedTree.Interval;
import augmentedTree.IntervalTree;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMRecord;

import java.util.*;

public class Gene implements Interval {
    private final int start;
    private final int end;
    private final String geneId;
    private final ArrayList<Transcript> transcriptList;
    private final String geneName;
    private final String bioType;
    private final String chr;
    private final char strand;
    private final ArrayList<Intron> introns = new ArrayList<>();
    private IntervalTree<Region> mappedAliBlocksTree = null;
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

    public void generateCDSIntrons() {
        for (Transcript transcript : transcriptList) {
            for (int i = 0; i < transcript.getCdsList().size() - 1; i++) {
                int intronStart = transcript.getCdsList().get(i).getStop() + 1;
                int intronEnd = transcript.getCdsList().get(i + 1).getStart() - 1;
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

    public void invertTranscriptsCds() {
        for (int i = 0; i < transcriptList.size(); i++) {
            Transcript currTranscript = transcriptList.get(i);
            currTranscript.reversCdsList();
            for (int j = 0; j < currTranscript.getCdsList().size(); j++) {
                currTranscript.getCdsList().get(j).setPos(j);
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
                                exonList.get(i).setTranscriptId(currTranscript.getTranscriptId());
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


    public ArrayList<Region> getSkippedCds() {
        // store events here
        ArrayList<Region> skippedCds = new ArrayList<>();
        boolean atLeastOneWT = false;
        // for every intron check every transcript
        for (Intron intron : introns) {
            int intronStart = intron.getStart();
            int intronEnd = intron.getEnd();

            for (Transcript currTranscript : transcriptList) {
                // get relevant HashMaps and check if currTranscript has cds starting or ending at i_s i_e
                HashMap<Integer, Region> cdsEnds = currTranscript.getCdsEnds();
                HashMap<Integer, Region> cdsStarts = currTranscript.getCdsStarts();

                boolean hasExonInFront = cdsEnds.containsKey(intronStart - 1);
                boolean hasExonBehind = cdsStarts.containsKey(intronEnd + 1);

                if (hasExonInFront && hasExonBehind) {
                    Region frontRegion = cdsEnds.get(intronStart - 1);
                    Region  behindRegion = cdsStarts.get(intronEnd + 1);

                    // get offset / look if there are exons inbetween
                    int offset = behindRegion.getPos() - frontRegion.getPos();

                    if (offset != 1) {
                        atLeastOneWT = true;
                        ArrayList<Region> cdsList = currTranscript.getCdsList();
                        for (int i = frontRegion.getPos(); i < behindRegion.getPos(); i++) {
                            if (i > frontRegion.getPos() && i < behindRegion.getPos()) {
                                cdsList.get(i).setTranscriptId(currTranscript.getTranscriptId());
                                skippedCds.add(cdsList.get(i));
                            }
                        }
                    }
                }
            }
        }
        if (atLeastOneWT) {
            return skippedCds;
        }
        return null;
    }
    public void addReadPairGap(Region fwLastBlock, Region rwFirstBlock) {
        ArrayList<Region> sorted = new ArrayList<>();
        sorted.add(fwLastBlock);
        sorted.add(rwFirstBlock);
        sorted.sort(Comparator.comparingInt((Region::getStart)));
        if (sorted.getFirst().getStop() < sorted.getLast().getStart()) {
            Region gap = new Region(sorted.getFirst().getId(), sorted.getFirst().getStop() + 1, sorted.getLast().getStart() - 1);
            // annotate gap
            gap.setTranscriptId(sorted.getFirst().getTranscriptId());
            gap.setType("GAP BETWEEN 2 READS");
            this.gappedAliBlocksTree.add(gap);
        }
    }

    public void addAlignedBlocks(TreeSet<Region> blocks, String type) {
        if (blocks.size() > 1) {
            int count = 1;
            int gapStart = 0;
            int gapEnd;
            for (Region r : blocks) {
                if (count % 2 != 0) {
                    if (gapStart == 0) {
                        gapStart = r.getStop() + 1;
                    } else {
                        gapEnd = r.getStart()-1;
                        Region gap = new Region(r.getId(), gapStart, gapEnd);
                        // annotate gap
                        gap.setTranscriptId(r.getTranscriptId());
                        gap.setType("GAP");
                        gappedAliBlocksTree.add(gap);
                        gapStart = r.getStop() + 1;

                    }
                    count++;
                } else {
                    gapEnd = r.getStart() - 1;
                    Region gap = new Region(r.getId(), gapStart, gapEnd);
                    // annotate gap
                    gap.setTranscriptId(r.getTranscriptId());
                    gap.setType("GAP");
                    gappedAliBlocksTree.add(gap);
                    gapStart = r.getStop() + 1;
                    count = 1;
                }
                Region mappedRead = new Region(r.getId(), r.getStart(), r.getStop());

                mappedRead.setTranscriptId(r.getTranscriptId());
                mappedRead.setType(type);
                this.mappedAliBlockSet.add(mappedRead);
            }
        }
        else { // there is only one alignment block and no gaps
            Region a = blocks.getLast();
            Region mappedRead = new Region(a.getId(), a.getStart(), a.getStop());
            mappedRead.setTranscriptId(a.getTranscriptId());
            mappedRead.setType(type);
            this.mappedAliBlockSet.add(mappedRead);
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
}

