package org.broadinstitute.quartetibdanalysis.tools;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.primitives.Ints;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.IntervalTree;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.MathUtils;
import scala.Tuple2;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

@CommandLineProgramProperties(
        summary = "Compute DD location based on IBD",
        oneLineSummary = "Compute DD location based on IBD",
        programGroup = StructuralVariantDiscoveryProgramGroup.class)
public class TargetDispersedDuplication extends GATKTool {

    @Argument(shortName = "ibdQuals", fullName = "ibdQuals", doc = "File of IBD quals", optional = true)
    private File ibdQuals = null;

    @Argument(shortName = "ibdQualsReduced", fullName = "ibdQualsReduced", doc = "File of IBD quals in reduced format", optional = true)
    private File ibdQualsReduced = null;

    @Argument(shortName = "annotatedCN", fullName = "annotatedCN", doc = "File of annotated copy numbers", optional = false)
    private File annotatedCN = null;

    @Argument(shortName = "outDir", fullName = "outDir", doc = "Output directory", optional = false)
    private File outDir = null;

    @Argument(shortName = "siblingIBD", fullName = "siblingIBD", doc = "Sibling IBD file", optional = false)
    private File siblingIBD = null;

    @Argument(shortName = "outputDataFile", fullName = "outputDataFile", doc = "Output reduced ibqq file", optional = true)
    private File outputDataFile = null;

    @Argument(shortName = "pairListFile", fullName = "pairListFile", doc = "List of pairs to limit results to", optional = true)
    private File pairListFile = null;

    @Argument(shortName = "cnvListFile", fullName = "cnvListFile", doc = "List of cnv ids to limit results to", optional = true)
    private File cnvListFile = null;

    @Argument(shortName = "excludedSamplesFile", fullName = "excludedSamplesFile", doc = "List samples to exclude from analysis", optional = true)
    private File excludedSamplesFile = null;

    @Argument(fullName = "lodScoreResolution", doc = "Min difference in LOD score at which to break intervals in output file ", optional = true)
    private double lodScoreResolution = 0.01;

    public static final double LOG_10_25_PCT = Math.log10(.25);
    public static final double LOG_10_50_PCT = Math.log10(.5);

    Map<String, Long> contigStarts = new HashMap<>();
    private Map<String, Tuple2<int[], Map<String, double[]>[]>> dataByChrom;
    private Map<String, Map<String, IntervalTree<Double>>> ibdOneIntervalTree;
    private Set<String> pairs;
    private Set<String> cnvIds;
    private Set<String> pairsWithIBDQData;
    private Set<String> excludedSamples;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    protected void onStartup() {

        super.onStartup();

        long pos = 0;
        for (SAMSequenceRecord sequence : getReferenceDictionary().getSequences()) {
            contigStarts.put(sequence.getSequenceName(), pos);
            pos = pos + sequence.getSequenceLength();
        }

        if (ibdQuals != null && ibdQualsReduced != null || (ibdQuals == null && ibdQualsReduced == null)) {
            throw new GATKException("need to specify one of ibdQuals, ibdQualsReduced");
        }

        try {

            pairsWithIBDQData = new HashSet<>();

            if (ibdQuals != null) {
                dataByChrom = readIBDData(ibdQuals);
            } else {
                dataByChrom = readReducedIBDData(ibdQualsReduced);
            }

            if (outputDataFile != null) {
                writeData(dataByChrom, outputDataFile);
            }
            ibdOneIntervalTree = readIDBOneRegions(siblingIBD);
            if (pairListFile != null) {
                pairs = readList(pairListFile);
            }

            if (cnvListFile != null) {
                cnvIds = readList(cnvListFile);
            }

            if (excludedSamplesFile != null) {
                excludedSamples = readList(excludedSamplesFile);
            } else {
                excludedSamples = new HashSet<>();
            }

        } catch (IOException e) {
            throw new GATKException(e.getMessage(), e);
        }

    }

    private Set<String> readList(final File listFile) throws IOException {
        Set<String> items = new HashSet<>();
        try(BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(listFile)))) {
            for (String line; (line = br.readLine()) != null; ) {
                items.add(line);
            }
        }
        return items;
    }

    @SuppressWarnings("unchecked")
    private Map<String, Tuple2<int[], Map<String, double[]>[]>> readReducedIBDData(final File ibdQualsReduced) throws IOException {
        Map<String, Tuple2<int[], Map<String, double[]>[]>> data = new HashMap<>();
        String prevChrom = "NA";
        List<Integer> positions = new ArrayList<>();
        List<Map<String, double[]>> mapsForChrom = new ArrayList<>();
        try(BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(ibdQualsReduced)), StandardCharsets.US_ASCII))) {
            for (String line; (line = br.readLine()) != null; ) {
                final String[] fields = line.split("\t");
                final String chrom = fields[0];
                if (!chrom.equals(prevChrom)) {
                    if (!"NA".equals(prevChrom)) {
                        System.err.println("Done processing " + prevChrom);
                        final int[] posArr = Ints.toArray(positions);
                        System.err.println("positions: " + posArr.length);
                        @SuppressWarnings("rawtypes")
                        final HashMap[] valArr = mapsForChrom.toArray(new HashMap[mapsForChrom.size()]);
                        data.put(prevChrom, new Tuple2<>(posArr, valArr));

                        positions = new ArrayList<>();
                        mapsForChrom = new ArrayList<>();
                    }
                    prevChrom = chrom;
                }
                final int pos = Integer.valueOf(fields[1]);
                positions.add(pos);
                final HashMap<String, double[]> vals = new HashMap<>(fields.length - 2);
                for (int i = 2; i < fields.length; i++) {
                    String[] mapFields = fields[i].split(":");
                    if (mapFields.length < 2) {
                        throw new GATKException("Bad line " + line + "; map fields: " + fields[i]);
                    }
                    String pair = mapFields[0];
                    double[] ibdqs = parseIbqs(mapFields[1]);
                    vals.put(pair, ibdqs);
                    pairsWithIBDQData.add(pair);
                }
                mapsForChrom.add(vals);

            }
        }

        System.err.println("Done processing " + prevChrom);

        final int[] posArr = Ints.toArray(positions);
        @SuppressWarnings("rawtypes")
        final HashMap[] valArr = mapsForChrom.toArray(new HashMap[mapsForChrom.size()]);
        System.err.println("position[5]: " + posArr[5]);
        System.err.println("valArr[5]: " + valArr[5]);

        data.put(prevChrom, new Tuple2<>(posArr, valArr));

        return data;
    }

    private double[] parseIbqs(final String mapField) {
        final double[] result = new double[3];
        final String[] ibqStrings = mapField.split(",");
        for (int i = 0; i < ibqStrings.length; i++) {
            String s = ibqStrings[i];
            result[i] = Double.valueOf(s);
        }
        return result;
    }

    private void writeData(final Map<String, Tuple2<int[], Map<String, double[]>[]>> dataByChrom, final File outputDataFile) throws IOException {
        final FileWriter fileWriter = new FileWriter(outputDataFile);

        for (final String contig : getReferenceDictionary().getSequences().stream().map(SAMSequenceRecord::getSequenceName).collect(Collectors.toList())) {
            if (dataByChrom.containsKey(contig)) {
                final Tuple2<int[], Map<String, double[]>[]> dataForChrom = dataByChrom.get(contig);
                final int[] breaks = dataForChrom._1();
                final Map<String, double[]>[] maps = dataForChrom._2();

                for (int i = 0; i < breaks.length; i++) {
                    int pos = breaks[i];
                    fileWriter.write(contig);
                    fileWriter.write("\t");
                    fileWriter.write(String.valueOf(pos));
                    for (String key : maps[i].keySet()) {
                        fileWriter.write("\t" + key + ":" + formatQuals(maps[i].get(key)));
                    }
                    fileWriter.write("\n");
                }
            }
        }
        fileWriter.close();
    }

    private String formatQuals(final double[] vals) {
        return vals[0] + "," + vals[1] + "," + vals[2];
    }

    private Map<String, Map<String, IntervalTree<Double>>> readIDBOneRegions(File siblingIBDFile) throws IOException {
        Map<String, Map<String, IntervalTree<Double>>> result = new HashMap<>();
        int linesRead = 0;
        try (BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(siblingIBDFile)), StandardCharsets.US_ASCII))) {
            for (String line; (line = br.readLine()) != null; ) {
                if (linesRead == 0) {
                    linesRead++;
                    continue;
                }
                final String[] fields = line.split("\t");
                final String pairName = fields[0];
                final String chrom = fields[3];
                final int start = Integer.parseInt(fields[4]);
                final int end = Integer.parseInt(fields[5]);
                final String state = fields[6];
                if ("ONE".equals(state)) {
                    // adding pseudocounts
                    final int patInd = Integer.parseInt(fields[7]) + 1;
                    final int matInd = Integer.parseInt(fields[8]) + 1;
                    final double patRatio = (double) patInd / (patInd + matInd);
                    final HashMap<String, IntervalTree<Double>> genomeIntervalTree = new HashMap<>();
                    final IntervalTree<Double> chromIntervalTree = new IntervalTree<>();
                    chromIntervalTree.put(start, end, patRatio);
                    if (! genomeIntervalTree.containsKey(chrom)) {
                        genomeIntervalTree.put(chrom, chromIntervalTree);
                    }
                    if (! result.containsKey(pairName)) {
                        result.put(pairName, genomeIntervalTree);
                    }

                }
            }
        }
        return result;
    }

    @SuppressWarnings("unchecked")
    private Map<String, Tuple2<int[], Map<String, double[]>[]>> readIBDData(File inputFile) throws IOException {
        final Map<String, Tuple2<int[],Map<String,double[]>[]>> dataByChrom = new HashMap<>(getReferenceDictionary().getSequences().size());

        final int chrom1Capacity = getReferenceDictionary().getSequence(0).getSequenceLength() / 200;
        List<Integer> positions = new ArrayList<>(chrom1Capacity);
        List<Map<String,double[]>> vals = new ArrayList<>(chrom1Capacity);
        //List<Tuple2<Long,>> data = new ArrayList<>(getReferenceDictionary().getSequence(0).getSequenceLength() / 100);
        Map<String, double[]> currentVals = new HashMap<>(800);
        int breaks = 0;
        long lastUpdate = 0;
        long lastMillis = System.currentTimeMillis();

        String prevChrom = "NA";
        try(BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputFile)), StandardCharsets.US_ASCII))) {
            for(String line; (line = br.readLine()) != null; ) {
                final String[] fields = line.split("\t");
                final String chrom = fields[0];
                if (! chrom.equals(prevChrom)) {
                    if (! "NA".equals(prevChrom)) {
                        System.err.println("Done processing " + prevChrom);
                        final int[] posArr = Ints.toArray(positions);
                        @SuppressWarnings("rawtypes")
                        final HashMap[] valArr = vals.toArray(new HashMap[vals.size()]);
                        System.err.println("Positions: " + posArr.length + ", vals: " + valArr.length);
                        dataByChrom.put(prevChrom, new Tuple2<>(posArr, valArr));
                        breaks = 0;

//                        for (int i = 0; i < posArr.length; i++) {
//                            System.out.print(prevChrom + "\t" + posArr[i]);
//                            final HashMap<String, double[]> pointVals = (HashMap<String, double[]>) valArr[i];
//                            for (String pair : pointVals.keySet()) {
//                                System.out.print("\t" + pair + "\t" + formatQuals(pointVals.get(pair)));
//                            }
//                            System.out.print("\n");
//                            if (i > 25) break;
//                        }
//
//                        for (int i = 0; i < posArr.length; i++) {
//                            final HashMap<String, double[]> pointVals = (HashMap<String, double[]>) valArr[i];
//                            if (pointVals.containsKey("SSC04928-SSC12955")) {
//                                System.out.println(prevChrom + "\t" + posArr[i] + "\t" + formatQuals(pointVals.get("SSC04928-SSC12955")));
//                            }
//                        }


                        vals = new ArrayList<>(getReferenceDictionary().getSequence(chrom).getSequenceLength() / 200);
                        positions = new ArrayList<>(getReferenceDictionary().getSequence(chrom).getSequenceLength() / 200);
                        currentVals = new HashMap<>(800);

                        lastUpdate = 0;

                        System.gc();
                        //Getting the runtime reference from system
                        Runtime runtime = Runtime.getRuntime();

                        int mb = 1024*1024;
                        System.out.println("##### Heap utilization statistics [MB] #####");

                        //Print used memory
                        System.out.println("Used Memory:"
                                + (runtime.totalMemory() - runtime.freeMemory()) / mb);

                        //Print free memory
                        System.out.println("Free Memory:"
                                + runtime.freeMemory() / mb);

                        //Print total available memory
                        System.out.println("Total Memory:" + runtime.totalMemory() / mb);

                        //Print Maximum available memory
                        System.out.println("Max Memory:" + runtime.maxMemory() / mb);
                    }
                }
                prevChrom = chrom;
                int pos = Integer.parseInt(fields[1]);
                String pair = fields[2];

                final String[] qualStrings = fields[4].split(",");
                if (qualStrings.length != 3) {
                    throw new GATKException("bad quals on line" + line);
                }
                final double[] quals = new double[3];
                quals[0] = Double.parseDouble(qualStrings[0]);
                quals[1] = Double.parseDouble(qualStrings[1]);
                quals[2] = Double.parseDouble(qualStrings[2]);

                //System.err.println(line);

                final double[] currentValForPair = currentVals.get(pair);
                if (currentValForPair == null || different(quals, currentValForPair)) {
                    if (breaks > 0 && positions.get(breaks - 1) == pos) {
                        //System.err.println("updating");
                        vals.get(breaks - 1).put(pair, quals);
                    } else {
                        final HashMap<String, double[]> locData = new HashMap<>(5);
                        //System.err.println("adding");
                        locData.put(pair, quals);
                        positions.add(pos);
                        vals.add(locData);
                        breaks = breaks + 1;
                    }
                    currentVals.put(pair, quals);
                }
                pairsWithIBDQData.add(pair);
                if (pos - lastUpdate > 1000000) {
                    final long newMillis = System.currentTimeMillis();
//                    System.err.println("position " + pos + "\ttime=" + (newMillis - lastMillis) + "\tbreaks="+breaks);
                    lastMillis = newMillis;
                    lastUpdate = pos;
                }
            }
        }

        final int[] posArr = Ints.toArray(positions);
        @SuppressWarnings("rawtypes")
        final HashMap[] valArr = vals.toArray(new HashMap[vals.size()]);

        dataByChrom.put(prevChrom, new Tuple2<>(posArr, valArr));

        return dataByChrom;
    }

    @VisibleForTesting
    private static boolean different(final double[] q1, final double[] q2) {
        //return (quals[0] != currentValForPair[0]) || (quals[1] != currentValForPair[1]) || (quals[2] != currentValForPair[2]);
        double q1Max = 0;
        double q1Min = Double.MAX_VALUE;
        int q1MaxIdx = 0;
        int q1MinIdx = 0;
        int zeroIdx = 0;
        for (int i = 0; i < q1.length; i++) {
            final double i1 = q1[i];
            if (i1 != 0 && i1 > q1Max) {
                q1Max = i1;
                q1MaxIdx = i;
            }
            if (i1 != 0 && i1 <= q1Min) {
                q1Min = i1;
                q1MinIdx = i;
            }
            if (i1 == 0) {
                zeroIdx = i;
            }
        }

        return q2[zeroIdx] != 0 || Math.abs((int) q1Min - (int) q2[q1MinIdx]) >= 5;
    }

    private long getLoc(final String chrom, final long pos) {
        return contigStarts.get(chrom) + pos;
    }

    @Override
    public void traverse() {

        String prevID = "NA";
        int linesRead = 0;
        List<QuartetCopyNumberEvent> informativeCopyNumberEvents = new ArrayList<>();
        int discordantCopyNumberEvents = 0;
        double expectedTrueDiscordantEvents = 0;
        Map<String, Integer> columnToField = new HashMap<>();
        try (BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(annotatedCN)), StandardCharsets.US_ASCII))) {
            for(String line; (line = br.readLine()) != null; ) {
                if (linesRead == 0) {
                    final String[] fields = line.split("\t");
                    for (int i = 0; i < fields.length; i++) {
                        columnToField.put(fields[i], i);
                    }
                    linesRead++;
                    continue;
                }
                final String[] fields = line.split("\t");

                final String sib1 = fields[columnToField.get("SIB1")];
                final String sib2 = fields[columnToField.get("SIB2")];
                final String father = fields[columnToField.get("PAT_ID")];
                final String mother = fields[columnToField.get("MAT_ID")];

                if (excludedSamples.contains(sib1) ||
                    excludedSamples.contains(sib2) ||
                    excludedSamples.contains(father) ||
                    excludedSamples.contains(mother)) {
                    continue;
                }

                QuartetCopyNumberEvent quartetCopyNumberEvent = new QuartetCopyNumberEvent(fields, columnToField);

                if (cnvIds != null && ! cnvIds.contains(quartetCopyNumberEvent.getID()) ) {
                    continue;
                }

                if (pairs != null && ! pairs.contains(quartetCopyNumberEvent.getPair()) ) {
                    continue;
                }

                if (! prevID.equals(quartetCopyNumberEvent.getID())) {

                    if (expectedTrueDiscordantEvents >= 2 && informativeCopyNumberEvents.size() >= 10) {
                        System.out.println("Processing " + prevID);
                        System.out.println("Discordant " + discordantCopyNumberEvents);
                        System.out.println("Expected true discordant " + expectedTrueDiscordantEvents);
                        System.out.println("Informative " + informativeCopyNumberEvents.size());
                        processInformativeEvents(informativeCopyNumberEvents,
                                discordantCopyNumberEvents,
                                expectedTrueDiscordantEvents,
                                lodScoreResolution);
                    }

                    informativeCopyNumberEvents = new ArrayList<>();
                    discordantCopyNumberEvents = 0;
                    expectedTrueDiscordantEvents = 0;
                    prevID = quartetCopyNumberEvent.getID();
                }

                if ( !quartetCopyNumberEvent.isError() && quartetCopyNumberEvent.allGenotyped() && !quartetCopyNumberEvent.isDeletion()
                        && quartetCopyNumberEvent.informative()) {
                    if (! pairsWithIBDQData.contains(quartetCopyNumberEvent.getPair())) {
                        System.err.println("No IBDQ data for potentially informative event " + quartetCopyNumberEvent.getID() + " / " + quartetCopyNumberEvent.getPair());
                    } else {
                        informativeCopyNumberEvents.add(quartetCopyNumberEvent);
                        if (!quartetCopyNumberEvent.isConsistent()) {
                            System.out.println("Discordant event pWrong: " + quartetCopyNumberEvent.getpWrong());
                            discordantCopyNumberEvents = discordantCopyNumberEvents + 1;
                            expectedTrueDiscordantEvents += 1 - quartetCopyNumberEvent.getpWrong();
                        }
                    }
                }

                linesRead++;
            }

            if (expectedTrueDiscordantEvents > 2 && informativeCopyNumberEvents.size() >= 10) {
                System.out.println("Processing " + prevID);
                System.out.println("Discordant " + discordantCopyNumberEvents);
                System.out.println("Informative " + informativeCopyNumberEvents.size());
                processInformativeEvents(informativeCopyNumberEvents,
                        discordantCopyNumberEvents,
                        expectedTrueDiscordantEvents,
                        lodScoreResolution);
            }

        } catch (FileNotFoundException e) {
            throw new GATKException("Could not open file", e);
        } catch (IOException e) {
            throw new GATKException("Could not read file", e);
        }
    }

    private void processInformativeEvents(final List<QuartetCopyNumberEvent> informativeCopyNumberEvents,
                                          final int discordantCopyNumberEvents,
                                          final double expectedTrueDiscordantEvents,
                                          final double lodScoreResolution) throws FileNotFoundException {

        Map<String, Tuple2<int[], double[]>> resultsByChrom = new HashMap<>();

        final QuartetCopyNumberEvent firstEvent = informativeCopyNumberEvents.get(0);
        final String cnvId = firstEvent.getID();

        final String cnvContig;
        final int cnvStart;


        final String[] cnvFields = cnvId.split("_");
        if (cnvId.startsWith("GS_SD_M2")) {
            cnvContig = cnvFields[3];
            cnvStart = Integer.parseInt(cnvFields[4]);
        } else {
            cnvContig = cnvFields[1];
            cnvStart = Integer.parseInt(cnvFields[2]);
        }

        boolean foundLoc =false;
        double cnvProb = Double.NEGATIVE_INFINITY;

        for (String contig : dataByChrom.keySet()) {
            Map<String, double[]> currentVals = new HashMap<>(800);

            System.err.println("processing contig " + contig);
            final Tuple2<int[], Map<String, double[]>[]> dataForChrom = dataByChrom.get(contig);
            final int[] breaks = dataForChrom._1();
            final Map<String, double[]>[] quals = dataForChrom._2();

            final double[] segmentScores = new double[breaks.length];

            int pos = 0;
            while (pos < breaks.length - 1) {
                segmentScores[pos] = Double.NEGATIVE_INFINITY;
                final Map<String, double[]> valsAtPos = quals[pos];
                currentVals.putAll(valsAtPos);

                double p = 0;

//                boolean debugSeg = "7".equals(contig) && breaks[pos] == 284829;
                boolean debugSeg = false;

                for (QuartetCopyNumberEvent cn : informativeCopyNumberEvents) {
                    double[] ibdq;
                    //double[] ibdqD;

                    Double patRatio;

                    if (debugSeg) System.err.println("pair: " + cn.getPair());
                    final String pair = cn.getPair();

                    if (currentVals.containsKey(pair)) {
                        ibdq = currentVals.get(pair);
                        //ibdqD = toExpDoubles(ibdq);
                        patRatio = getPatRatio(contig, breaks, pos, pair);
                    } else {
                        ibdq = getMissingIBDQValues(contig, breaks, quals, pos, debugSeg, pair);
                        patRatio = null;
                    }



                    double cnp = cn.pHere(ibdq, patRatio, debugSeg);

                    if (debugSeg) System.err.println("cnp: " + cnp);
                    p = p + cnp;
                    if (debugSeg) System.err.println("p: " + p);

                }

                if (! foundLoc && contig.equals(cnvContig) && breaks[pos] < cnvStart && breaks[pos+1] >= cnvStart) {
                    foundLoc = true;
                    cnvProb = p;
                }

                segmentScores[pos] = p;
                pos = pos + 1;
            }
            resultsByChrom.put(contig, new Tuple2<>(breaks, segmentScores));
        }

        double maxLod = Double.NEGATIVE_INFINITY;
        String maxLodChr = "";
        int maxLodStart = -1;
        int maxLodEnd = -1;
        for (String contig : dataByChrom.keySet()) {
            final Tuple2<int[], double[]> resultsForChrom = resultsByChrom.get(contig);
            final double[] scores = resultsForChrom._2();
            for (int i = 0; i < scores.length - 1; i++) {
                scores[i] = scores[i] - cnvProb;
                if (scores[i] > maxLod) {
                    maxLod = scores[i];
                    maxLodChr = contig;
                    maxLodStart = resultsForChrom._1()[i];
                    maxLodEnd = resultsForChrom._1()[i + 1];
                }
            }
        }


        try (PrintStream printStream = new PrintStream(new BufferedOutputStream(new GZIPOutputStream(new FileOutputStream(new File(outDir + "/" + cnvId + ".bedgraph.gz")))))) {

            for (final String contig : getReferenceDictionary().getSequences().stream().map(SAMSequenceRecord::getSequenceName).collect(Collectors.toList())) {
                if (! resultsByChrom.containsKey(contig)) {
                    continue;
                }
                final Tuple2<int[], double[]> resultsForChrom = resultsByChrom.get(contig);
                final int[] breaks = resultsForChrom._1();
                final double[] segmentScores = resultsForChrom._2();
                double prevScore = 1;
                int prevStart = -1;
                for (int i = 0; i < breaks.length - 1; i++) {
                    final int start = breaks[i];
                    final double score = segmentScores[i];
                    if (Math.abs(score - prevScore) > lodScoreResolution) {
                        if (prevStart != -1) {
                            printStream.println(contig + "\t" + prevStart + "\t" + start + "\t" + prevScore);
                        }
                        prevStart = start;
                        prevScore = score;
                    }
                }
            }
        } catch (IOException e) {
            throw new GATKException(e.getMessage(), e);
        }

        try (PrintStream printStream = new PrintStream(new BufferedOutputStream(new FileOutputStream(new File(outDir + "/" + cnvId + ".metadata"))))) {
            printStream.println("ID\tDISCORDANT\tEXPECTED_DISCORDANT\tINFORMATIVE\tORIGINAL_PROB\tMAX_LOD\tMAX_LOD_CHR\tMAX_LOD_START\tMAX_LOD_END");
            printStream.println(cnvId + "\t" + discordantCopyNumberEvents + "\t" + expectedTrueDiscordantEvents + "\t" + informativeCopyNumberEvents.size() + "\t" + cnvProb + "\t" + maxLod + "\t" + maxLodChr + "\t" + maxLodStart + "\t" + maxLodEnd);
        } catch (IOException e) {
            throw new GATKException(e.getMessage(), e);
        }

        calculateConfidenceIntervals(resultsByChrom, cnvId, cnvProb);
    }

    private void calculateConfidenceIntervals(final Map<String, Tuple2<int[], double[]>> resultsByChrom, final String cnvId, final double cnvProb) {


        final int numSegments = resultsByChrom.values().stream().mapToInt(pv -> pv._2().length).sum();

        final List<MapSegment> segments = new ArrayList<>(numSegments);

        for (String contig : resultsByChrom.keySet()) {
            final Tuple2<int[], double[]> resultsForContig = resultsByChrom.get(contig);
            final int length = resultsForContig._1().length;
            for (int i = 0; i < length - 1; i++) {
                segments.add(new MapSegment(contig, resultsForContig._1()[i], resultsForContig._1[i+1], resultsForContig._2[i] + cnvProb));
            }
            segments.add(new MapSegment(contig, resultsForContig._1()[length - 2], resultsForContig._1()[length - 1], resultsForContig._2()[length - 1] + cnvProb));
        }

        segments.sort(Comparator.comparing(m -> m.likelihood));

        final double maxLikelihood =
                segments
                        .stream()
                        .mapToDouble(s -> s.likelihood)
                        .max().orElse(0);

        System.err.println("max likelihood: " + maxLikelihood);

        final double genomeSize =
                segments
                        .stream()
                        .mapToDouble(s -> s.end - s.start)
                        .max().orElse(0);

        System.err.println("genome size: " + genomeSize);

        double totalLikelihood =
                    maxLikelihood +
                    StrictMath.log10(
                            segments
                                    .stream()
                                    .mapToDouble(s ->  StrictMath.pow(10, s.likelihood - maxLikelihood) * (s.end - s.start) / genomeSize)
                                    .sum());

        System.err.println("total likelihood: " + totalLikelihood);

        final List<MapSegment> confidenceInterval = new ArrayList<>();
        double confidenceProb = 0;

        int j = segments.size();

        final Set<Double> writtenCIs = new HashSet<>();

        final Comparator<MapSegment> comparator = new MapSegmentComparator(getReferenceDictionary());

        while (confidenceProb < 0.99) {
            final MapSegment nextSegment = segments.get(j - 1);
            final double nextProb = StrictMath.pow(10, nextSegment.likelihood - totalLikelihood) * (nextSegment.end - nextSegment.start) / genomeSize;
            confidenceProb = confidenceProb + nextProb;
            System.err.println("adding segment to CI " + nextSegment.contig + ":" + nextSegment.start + "-" + nextSegment.end + "\t" + nextSegment.likelihood + "\t" + nextProb + " (" + confidenceProb + ")");
            confidenceInterval.add(nextSegment);
            if (confidenceProb > 0.95 && ! writtenCIs.contains(0.95)) {
                confidenceInterval.sort(comparator);
                try (PrintStream printStream95 = new PrintStream(new BufferedOutputStream(new FileOutputStream(new File(outDir + "/" + cnvId + ".95ci.bed"))))) {
                    confidenceInterval.forEach(s -> printStream95.println(s.contig + "\t" + s.start + "\t" + s.end + "\t" + (StrictMath.pow(10, s.likelihood - totalLikelihood) * (nextSegment.end - nextSegment.start) / genomeSize)));
                } catch (IOException e) {
                    throw new GATKException(e.getMessage(), e);
                }
                writtenCIs.add(0.95);
            }
            if (confidenceProb > 0.99 && ! writtenCIs.contains(0.99)) {
                confidenceInterval.sort(comparator);
                try (PrintStream printStream99 = new PrintStream(new BufferedOutputStream(new FileOutputStream(new File(outDir + "/" + cnvId + ".99ci.bed"))))) {
                    confidenceInterval.forEach(s -> printStream99.println(s.contig + "\t" + s.start + "\t" + s.end + "\t" + (StrictMath.pow(10, s.likelihood - totalLikelihood) * (nextSegment.end - nextSegment.start) / genomeSize)));
                } catch (IOException e) {
                    throw new GATKException(e.getMessage(), e);
                }
                writtenCIs.add(0.99);
            }

            j--;
        }

    }

    static class MapSegment {
        String contig;
        int start;
        int end;
        double likelihood;

        public MapSegment(final String contig, final int start, final int end, final double likelihood) {
            this.contig = contig;
            this.start = start;
            this.end = end;
            this.likelihood = likelihood;
        }

        public String getContig() {
            return contig;
        }

        public int getStart() {
            return start;
        }

        public int getEnd() {
            return end;
        }

        public double getLikelihood() {
            return likelihood;
        }

    }

    static class MapSegmentComparator implements Comparator<MapSegment> {

        private final SAMSequenceDictionary referenceDictonary;

        public MapSegmentComparator(final SAMSequenceDictionary referenceDictionary) {
            this.referenceDictonary = referenceDictionary;
        }

        @Override
        public int compare(final MapSegment o1, final MapSegment o2) {
            int contig1 = referenceDictonary.getSequenceIndex(o1.contig);
            int contig2 = referenceDictonary.getSequenceIndex(o2.contig);
            if (contig1 < contig2) return -1;
            if (contig1 > contig2) return 1;
            return Integer.compare(o1.start, o2.start);
        }
    }

    private double[] getMissingIBDQValues(final String contig, final int[] breaks, final Map<String, double[]>[] quals, final int pos, final boolean debugSeg, final String pair) {
        //System.err.println("Segment " + contig + ":" + breaks[pos] + ": no IBD data available for pair " + pair + ", looking ahead");
        int pos2 = pos + 1;
        double[] ibdq = new double[] {LOG_10_25_PCT, LOG_10_50_PCT, LOG_10_25_PCT};
        while (pos2 < breaks.length) {
            if (quals[pos2].containsKey(pair)) {
                ibdq = quals[pos2].get(pair);
                //ibdqD = toExpDoubles(ibdq);
                break;
            }
            pos2++;
        }
        return ibdq;
    }

    private Double getPatRatio(final String contig, final int[] breaks, final int pos, final String pair) {
        Double patRatio = null;
        if (ibdOneIntervalTree.containsKey(pair)) {
            final Map<String, IntervalTree<Double>> genomeIbd1Map = ibdOneIntervalTree.get(pair);
            if (genomeIbd1Map.containsKey(contig)) {
                final IntervalTree.Node<Double> ibd1Node = genomeIbd1Map.get(contig).find(breaks[pos], breaks[pos + 1]);
                if (ibd1Node != null) {
                    patRatio = ibd1Node.getValue();
                }
            }
        }
        return patRatio;
    }


    private double[] toExpDoubles(final double[] ibdq) {
        final double[] result = new double[3];
        int zeroCount = 0;
        int zeroIdx1=-1;
        int zeroIdx2=-1;
        double total = 0;
        for (int i = 0; i < 3; i++) {
             if (ibdq[i] == 0) {
                 zeroCount = zeroCount + 1;
                 if (zeroIdx1 == -1) {
                     zeroIdx1 = i;
                 } else if (zeroIdx2 == -1) {
                     zeroIdx2 = i;
                 } else {
                     // three zeroes
                     result[0] = 1.0 / 3;
                     result[1] = 1.0 / 3;
                     result[2] = 1.0 / 3;
                 }
             } else {
                 result[i] = Math.pow(10, ibdq[i]);
                 total = total + result[i];
             }
        }
        if (zeroCount == 1) {
            result[zeroIdx1] = 1.0 - total;
        } else if (zeroCount == 2) {
            result[zeroIdx1] = (1.0 - total) / 2;
            result[zeroIdx2] = (1.0 - total) / 2;
        }
        return result;
    }

    static class QuartetCopyNumberEvent {

        private final String id;
        private final String pair;
        private final int sib1cn;
        private final double sib1cnq;
        private final int sib2cn;
        private final double sib2cnq;
        private final int patcn;
        private final double patcnq;
        private final int matcn;
        private final double matcnq;
        private final boolean constIbd0;
        private final boolean constIbd1;
        private final boolean constIbd1p;
        private final boolean constIbd1m;
        private final boolean constIbd2;
        private final boolean stateConst;
        private final String errorCond;
        private final double pWrong;

        private final double totalProbIbd0;
        private final double totalProbIbd1;
        private final double totalProbIbd1p;
        private final double totalProbIbd1m;
        private final double totalProbIbd2;

        private final boolean valid;

        public QuartetCopyNumberEvent(final String[] fields, final Map<String, Integer> columnToField) {

            this.id = fields[columnToField.get("ID")];
            this.pair = fields[columnToField.get("SIB_PAIR")];

            this.constIbd0 = Boolean.parseBoolean(fields[columnToField.get("CONST_IBD0")]);
            this.constIbd1 = Boolean.parseBoolean(fields[columnToField.get("CONST_IBD1")]);
            this.constIbd1p = Boolean.parseBoolean(fields[columnToField.get("CONST_IBD1P")]);
            this.constIbd1m = Boolean.parseBoolean(fields[columnToField.get("CONST_IBD1M")]);
            this.constIbd2 = Boolean.parseBoolean(fields[columnToField.get("CONST_IBD2")]);
            this.stateConst = Boolean.parseBoolean(fields[columnToField.get("STATE_CONSISTENT")]);
            this.errorCond = fields[columnToField.get("ERROR_COND")];

            this.sib1cn = Integer.parseInt(fields[columnToField.get("SIB1_CN")]);
            this.sib1cnq = "NA".equals(fields[columnToField.get("SIB1_CNQ")]) ? 0 : Double.parseDouble(fields[columnToField.get("SIB1_CNQ")]);
            this.sib2cn = Integer.parseInt(fields[columnToField.get("SIB2_CN")]);
            this.sib2cnq = "NA".equals(fields[columnToField.get("SIB2_CNQ")]) ? 0 :Double.parseDouble(fields[columnToField.get("SIB2_CNQ")]);
            this.patcn = Integer.parseInt(fields[columnToField.get("PAT_CN")]);
            this.patcnq = "NA".equals(fields[columnToField.get("PAT_CNQ")]) ? 0 :Double.parseDouble(fields[columnToField.get("PAT_CNQ")]);
            this.matcn = Integer.parseInt(fields[columnToField.get("MAT_CN")]);
            this.matcnq = "NA".equals(fields[columnToField.get("MAT_CNQ")]) ? 0 :Double.parseDouble(fields[columnToField.get("MAT_CNQ")]);

            this.totalProbIbd0 = "NA".equals(fields[columnToField.get("TOTAL_PROB_IBD0")]) ? 0 :Double.parseDouble(fields[columnToField.get("TOTAL_PROB_IBD0")]);
            this.totalProbIbd1 = "NA".equals(fields[columnToField.get("TOTAL_PROB_IBD1")]) ? 0 :Double.parseDouble(fields[columnToField.get("TOTAL_PROB_IBD1")]);
            this.totalProbIbd1p = "NA".equals(fields[columnToField.get("TOTAL_PROB_IBD1P")]) ? 0 :Double.parseDouble(fields[columnToField.get("TOTAL_PROB_IBD1P")]);
            this.totalProbIbd1m = "NA".equals(fields[columnToField.get("TOTAL_PROB_IBD1M")]) ? 0 :Double.parseDouble(fields[columnToField.get("TOTAL_PROB_IBD1M")]);
            this.totalProbIbd2 = "NA".equals(fields[columnToField.get("TOTAL_PROB_IBD2")]) ? 0 :Double.parseDouble(fields[columnToField.get("TOTAL_PROB_IBD2")]);

            pWrong = calcPWrong(sib1cnq, sib2cnq, patcnq, matcnq);

            this.valid = true;

        }

        static double calcPWrong(final double sib1cnq, final double sib2cnq, final double patcnq, final double matcnq) {
            return 1 - (1 - Math.pow(10, sib1cnq / -10)) * (1 - Math.pow(10, sib2cnq / -10)) * (1 - Math.pow(10, patcnq / -10)) * (1 - Math.pow(10, matcnq / -10));
        }

        public boolean isError() {
            return !"NA".equals(errorCond);
        }

        public boolean allGenotyped() {
            return sib1cn != -1 && sib2cn != -1 && patcn != -1 && matcn != -1;
        }

        public boolean informative() {
            return (constIbd0 || constIbd1 || constIbd2) && (!constIbd0 || !constIbd1 || !constIbd2);
        }

        public boolean isConsistent() {
            return stateConst;
        }

        public double pHere(double[] ibdQ, Double patRatio, final boolean debugSeg) {
            double total = Double.NEGATIVE_INFINITY;

            // totalProbIbd0 = p(const|ibd0)
            // total = sum(states) p(const|state)p(state) = p(const)

            // System.err.println("sum ibdq: "+ (Math.pow(10, ibdQ[0]) + Math.pow(10, ibdQ[1]) + Math.pow(10, ibdQ[2])));

            total = log10SumLog10Capped(new double[] {total,
                    totalProbIbd0 + ibdQ[0],
                    getIBD1Total(ibdQ, patRatio, totalProbIbd1, totalProbIbd1p, totalProbIbd1m, constIbd1p, constIbd1m),
                    totalProbIbd2 + ibdQ[2]
            });

            return total;
        }

        private double getIBD1Total(final double[] ibdQ, final Double patRatio, final double totalProbIbd1, final double totalProbIbd1p, final double totalProbIbd1m, final boolean constIbd1p, final boolean constIbd1m) {
            if (patRatio == null) {
                return totalProbIbd1 + ibdQ[1];
            } else {
                return log10SumLog10Capped(new double[] {
                        totalProbIbd1p +  (Math.log10(patRatio) + ibdQ[1]),
                        totalProbIbd1m + (Math.log10(1 - patRatio) + ibdQ[1])
                });
            }
        }

        private double log10SumLog10Capped(final double[] log10Values) {
            return Math.min(0.0, MathUtils.log10SumLog10(log10Values));
        }

        public String getPair() {
            return pair;
        }

        public String getID() {
            return id;
        }

        public boolean isDeletion() {
            return sib1cn <= 2 && sib2cn <= 2 && patcn <= 2 && matcn <= 2;
        }

        public double getpWrong() {
            return pWrong;
        }
    }

}
