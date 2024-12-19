package org.src;

import augmentedTree.IntervalTree;
import net.sf.samtools.*;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class BamFeatures {

    private final SAMFileReader samReader;
    private final Genome genome;

    private final HashSet<Gene> mappedGenes = new HashSet();

    public BamFeatures(String pathToBAM, String pathToGTF) throws IOException {
        this.genome = new Genome();
        this.genome.readGTF(pathToGTF);
        this.samReader = new SAMFileReader(new File(pathToBAM), false);
        this.samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
    }

    public void processBAM(String outPath) throws IOException {
        Iterator<SAMRecord> it = samReader.iterator();
        HashMap<String, SAMRecord> seenEntries = new HashMap<>();
        File outFile = new File(outPath);
        File parentDir = outFile.getParentFile();
        parentDir.mkdirs();
        String currentChr = null;

        while (it.hasNext()) {
            SAMRecord current = it.next();

            if (currentChr == null) {
                currentChr = current.getReferenceName();
            } else if (!currentChr.equals(current.getReferenceName())) {
                // clear seen
                seenEntries.clear();
                // clear chr in intervaltree as well
                genome.getIntervalTreeMap().get(currentChr).clear();
                // update currChr
                currentChr = current.getReferenceName();
            }

            if (!flagCheck(current)) {
                continue;
            }

            // track entries
            if (!seenEntries.containsKey(current.getReadName())) {
                seenEntries.put(current.getReadName(), current);
                continue;
            }

            // at this point we already have the read pair
            // get mate
            SAMRecord mate = seenEntries.get(current.getReadName());
            // init readpair with its correct strandness
            Boolean frstrand = null;
//            ReadPair pair = determineReadPair(mate, current, mate.getMateNegativeStrandFlag());
            ReadPair pair;
            if (mate.getFirstOfPairFlag()) {
                pair = new ReadPair(mate, current, !mate.getMateNegativeStrandFlag());
            } else {
                pair = new ReadPair(current, mate, !current.getMateNegativeStrandFlag());
            }

            int cgenes = pair.getcgenes(genome);

            if (cgenes == 0) {
                int igenes = pair.getigenes(genome);
                if (igenes > 0) {
                    continue;
                }
            }

            // check for split inconsistency
            int nsplit = pair.getNsplit();
            if (nsplit == -1) {
                continue;
            }
            // update count based on annotation
            ArrayList<Gene> transcriptomicGenes = pair.getTranscriptomicGenes();
            // skip if empty
            if (transcriptomicGenes.isEmpty()) {
                continue;
            }
            mappedGenes.addAll(transcriptomicGenes);
        }
    }
    public void getPctSplicedCunts() {
        for (Gene g : mappedGenes) {
            if (g.getStrand() == '-') {
                g.invertTranscripts();
            }
            g.generateIntrons();

//            if (g.getGeneId().equals("ENSG00000158109.10")) {
//                System.out.println();
//            }
//            if (g.getGeneId().equals("ENSG00000142634.8")) {
//                System.out.println();
//            }
            if (g.getGeneId().equals("ENSG00000132881.7")) {
                System.out.println();
            }
//            if (g.getGeneId().equals("ENSG00000158195.6")) {
//                System.out.println();
//            }
//            if (g.getGeneId().equals("ENSG00000198483.8")) {
//                System.out.println();
//            }
//            if (g.getGeneId().equals("ENSG00000196407.7")) {
//                System.out.println();
//            }
//            if (g.getGeneId().equals("ENSG00000158195.6")) {
//                System.out.println();
//            }

            ArrayList<Exon> skippedExons = g.getSkippedExons();

            if (skippedExons == null) {
                continue;
            }

            IntervalTree<Region> mappedAliReadsTree = g.getMappedAliBlocksTree();
            IntervalTree<Region> gappedAliReadsTree = g.getGappedAliBlocksTree();

            for (Exon skippedExon : skippedExons) {
                ArrayList<Region> inclusionReads = mappedAliReadsTree.getIntervalsSpannedBy(skippedExon.getStart(), skippedExon.getStop(), new ArrayList<>());
                HashSet<Region> exclusionReads = gappedAliReadsTree.getIntervalsSpanning(skippedExon.getStart(), skippedExon.getStop(), new HashSet<>());

                // skipp unmapped exons
                if(inclusionReads.isEmpty()) {
                    continue;
                }
                int total = inclusionReads.size() + exclusionReads.size();
                System.out.println(g.getGeneId() + " " + skippedExon.getStart() + "-" + skippedExon.getStop()+ " " + inclusionReads.size() + " " + exclusionReads.size() + " " + total);
            }
        }
    }
    public boolean flagCheck(SAMRecord record) {
        // ignore based on flags
        boolean isPrimary = !record.getNotPrimaryAlignmentFlag();
        boolean isMateMapped = !record.getMateUnmappedFlag();
        boolean isMapped = !record.getReadUnmappedFlag();
        boolean sameChr = record.getReferenceName().equals(record.getMateReferenceName());
        boolean oppStrand = record.getReadNegativeStrandFlag() != record.getMateNegativeStrandFlag();
        boolean paired = record.getReadPairedFlag();
        return isPrimary && isMapped && isMateMapped && sameChr && oppStrand && paired;
    }

    public ReadPair determineReadPair(SAMRecord mate, SAMRecord current, Boolean frstrand) {
        // bases on frstrand, getFirstOfPair and getNegativeStrandFlag determine
        // correct readpair configuration
        if (frstrand == null) {
            if (mate.getFirstOfPairFlag()) {
                return new ReadPair(mate, current, null);
            } else {
                return new ReadPair(current, mate, null);
            }
        }
        // fr +
        else if (frstrand) {
            // mate first
            if (mate.getFirstOfPairFlag()) {
                // curr -
                if (current.getReadNegativeStrandFlag()) {
                    return new ReadPair(current, mate, false);
                }
                // curr +
                else {
                    return new ReadPair(current, mate, true);
                }
            } else {
                // mate -
                if (mate.getReadNegativeStrandFlag()) {
                    return new ReadPair(mate, current, false);
                }
                // mate +
                else {
                    return new ReadPair(mate, current, true);
                }
            }
        }
        // fr -
        else {
            if (mate.getFirstOfPairFlag()) {
                // mate -
                if (mate.getReadNegativeStrandFlag()) {
                    return new ReadPair(mate, current, false);
                }
                // mate +
                else {
                    return new ReadPair(mate, current, true);
                }
            } else {
                // curr -
                if (current.getReadNegativeStrandFlag()) {
                    return new ReadPair(current, mate, false);
                }
                // curr +
                else {
                    return new ReadPair(current, mate, true);
                }
            }
        }
    }
}
