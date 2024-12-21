package org.src;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;
import java.io.IOException;

public class Main {
    public static void main(String[] args) throws IOException {
        ArgumentParser parser = ArgumentParsers.newFor("PsiCounter").build().defaultHelp(true).description("Usage:\n\t-gtf <path-to-gtf>\n\t-o <path-to-out.psi>\n\t-bam <path-to-bam>");
        try {
            parser.addArgument("-gtf").required(true).help("Path to Gene Transfer Format File.");
            parser.addArgument("-bam").required(true).help("Path to Bam File.");
            parser.addArgument("-o").required(true).help("Specify Output File Name.");

            Namespace ns = parser.parseArgs(args);
            String gtfPath = ns.getString("gtf");
            String out = ns.getString("o");
            String bamPath = ns.getString("bam");

            BamFeatures bam = new BamFeatures(bamPath, gtfPath);
            bam.processBAM();
            bam.getPctSplicedCounts(out);
        } catch (ArgumentParserException e) {
            parser.printHelp();
        }
    }
}