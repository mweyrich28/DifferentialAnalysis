package org.src;

import augmentedTree.IntervalTree;
import com.sun.source.tree.Tree;
import net.sf.samtools.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class BamFeatures {

    private final SAMFileReader samReader;
    private final Genome genome;

    private final IntervalTree<Region> readBlocks = new IntervalTree<>();
    private final IntervalTree<Region> readSpaces = new IntervalTree<>();
    private final HashSet<String> mappedTranscripts = new HashSet<>();

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
        BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(outPath));
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
            ReadPair pair = determineReadPair(mate, current, frstrand);

            // append read id to sb
            StringBuilder sb = new StringBuilder(current.getReadName());
            int cgenes = pair.getcgenes(genome);
            int gdist = 0;

            if (cgenes == 0) {
                int igenes = pair.getigenes(genome);
                if (igenes > 0) {
                    continue;
                }
            }

            int nsplit = pair.getNsplit();
            if (nsplit == -1) {
                continue;
            }
            // update count based on annotation
            String annotation = pair.annotateRegion();
            cgenes = pair.getgCount();

        }

        for (Gene g : genome.getGenes()) {
            if (g.getStrand() == '-') {
                g.invertTranscripts();
            }
            g.generateIntrons();

//            if (g.getGeneId().equals("ENSG00000158109.10")) {
//                System.out.println();
//            }

            ArrayList<Exon> skippedExons = g.getSkippedExons();

            if (skippedExons == null) {
                continue;
            }
            for (Exon skippedExon : skippedExons) {
                System.out.println(g.getGeneId() + " " + skippedExon.getStart() + "-" + skippedExon.getStop());
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
