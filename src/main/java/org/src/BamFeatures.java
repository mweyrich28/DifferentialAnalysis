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

    public BamFeatures(String pathToBAM, String pathToGTF) throws IOException {
        this.genome = new Genome();
        this.genome.readGTF(pathToGTF);
        this.samReader = new SAMFileReader(new File(pathToBAM), false);
        this.samReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
    }

    public void processBAM(String outPath) throws IOException {
        Iterator<SAMRecord> it = samReader.iterator();
        File outFile = new File(outPath);
        File parentDir = outFile.getParentFile();
        parentDir.mkdirs();
        BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(outPath));
        while (it.hasNext()) {
            SAMRecord current = it.next();

            if (!flagCheck(current)) {
                continue;
            }

            // add reads to intervalltree
        }
        for (Gene g : genome.getGenes()) {
            if (g.getStrand() == '-') {
                g.invertTranscripts();
            }
            g.generateIntrons();

            if (g.getGeneId().equals("ENSG00000158109.10")) {
                System.out.println();
            }
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
}
