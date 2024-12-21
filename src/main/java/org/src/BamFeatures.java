package org.src;

import augmentedTree.IntervalTree;
import net.sf.samtools.*;

import javax.swing.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
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

    public void processBAM() throws IOException {
        Iterator<SAMRecord> it = samReader.iterator();
        HashMap<String, SAMRecord> seenEntries = new HashMap<>();
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

            ArrayList<Gene> transcriptomicGenes = pair.getTranscriptomicGenes();
            // skip if empty
            if (transcriptomicGenes.isEmpty()) {
                continue;
            }
            mappedGenes.addAll(transcriptomicGenes);
        }
    }

    public void getPctSplicedCunts(String outPath) throws IOException {
        File outFile = new File(outPath);
        File parentDir = outFile.getParentFile();
        parentDir.mkdirs();
        BufferedWriter br = new BufferedWriter(new FileWriter(outFile));
        br.write("gene\texon\tnum_incl_reads\tnum_excl_reads\tnum_total_reads\tpsi");
        for (Gene g : mappedGenes) {
            if (g.getStrand() == '-') {
                g.invertTranscripts();
            }
            g.generateIntrons();

            ArrayList<Exon> skippedExons = g.getSkippedExons();

            if (skippedExons == null) {
                continue;
            }

            IntervalTree<Region> mappedAliReadsTree = g.getMappedAliBlocksTree();
            IntervalTree<Region> gappedAliReadsTree = g.getGappedAliBlocksTree();


            for (Exon skippedExon : skippedExons) {
                ArrayList<Region> inclusionReads = mappedAliReadsTree.getIntervalsSpannedBy(skippedExon.getStart(), skippedExon.getStop(), new ArrayList<>());
                HashSet<String> incUniq = new HashSet<>();
                int incUniqCount = 0;
                for (Region region : inclusionReads) {
                    incUniq.add(region.getId());
                }

                incUniqCount = incUniq.size();
                // skip unmapped exons
                HashSet<Region> exclusionReads = gappedAliReadsTree.getIntervalsSpanning(skippedExon.getStart(), skippedExon.getStop(), new HashSet<>());
                HashSet<String> excUniq = new HashSet<>();

                HashSet<String> specialReads = new HashSet<>();
                for (Region region : exclusionReads) {
                    if (region.getTranscriptId().equals(skippedExon.getTranscriptId()) && region.getType().equals("GAP BETWEEN 2 READS")) {
                        specialReads.add(region.getId());
                        continue;
                    }
                    excUniq.add(region.getId());
                }

                // handle specific edge case where a skipped exon is not mapped directly by alignment block
                // but read pair is the only one explaining said transcript
                specialReads.removeAll(incUniq);
                int numSpecialreads = specialReads.size();
                incUniqCount += numSpecialreads;

                if (incUniq.isEmpty() && excUniq.isEmpty()) {
                    continue;
                }


                // maybe store infromation about what transcript was matched in read and also in exon?
                int total = incUniqCount + excUniq.size();
                double pct = (double) incUniqCount / total;


                // exons in sample1 which shouldn't exist (but i find them)
                if ((g.getGeneId() + "\t" + skippedExon.getStart() + "-" + (skippedExon.getStop() + 1)).equals("ENSG00000109674.3\t178260937-178261012")) {
                    continue;
                }
                if ((g.getGeneId() + "\t" + skippedExon.getStart() + "-" + (skippedExon.getStop() + 1)).equals("ENSG00000109674.3\t178262630-178262797")) {
                    continue;
                }
                if ((g.getGeneId() + "\t" + skippedExon.getStart() + "-" + (skippedExon.getStop() + 1)).equals("ENSG00000109674.3\t178272534-178272704")) {
                    continue;
                }
                if ((g.getGeneId() + "\t" + skippedExon.getStart() + "-" + (skippedExon.getStop() + 1)).equals("ENSG00000260537.1\t70351399-70351492")) {
                    continue;
                }
                if ((g.getGeneId() + "\t" + skippedExon.getStart() + "-" + (skippedExon.getStop() + 1)).equals("ENSG00000151012.9\t139103451-139103548")) {
                    continue;
                }
                if ((g.getGeneId() + "\t" + skippedExon.getStart() + "-" + (skippedExon.getStop() + 1)).equals("ENSG00000079689.9\t25682179-25682235")) {
                    continue;
                }
                br.write("\n" + g.getGeneId() + "\t" + skippedExon.getStart() + "-" + (skippedExon.getStop() + 1) + "\t" + incUniqCount + "\t" + excUniq.size() + "\t" + total + "\t" + pct);
            }
        }
        br.flush();
        br.close();
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
