package org.broadinstitute.quartetibdanalysis.tools;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;

import java.io.File;
import java.util.Collections;

@CommandLineProgramProperties(
        summary = "Annotate a copy number VCF with IBD state",
        oneLineSummary = "Annotate a copy number VCF with IBD state",
        programGroup = QuartetIBDAnalysisProgramGroup.class
)
public class AnnotateCopyNumberVariantsWIthIBDState extends VariantWalker {

    @Argument(doc = "uri for the output file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = true)
    public File out;

    @Argument(fullName = "filterGenotypes", doc = "True to ignore genotypes that have been filtered (default true)")
    private Boolean filterGenoypes = null;

    @Argument(fullName = "genotypeQualityThreshold", doc = "Ignore genotypes below this genotype quality GQ value (default no threshold)")
    private Double genotypeQualityThreshold = null;

    @Argument(fullName = "siblingIBDFile", doc = "input Sibling IBD table")
    private File siblingIBDFile = null;

    @Argument(fullName = "minCallRate", doc = "lower threshold for CNV call rate to be considered for annotation", optional = true)
    private Float minCallRate = null;

    @Argument(fullName = "endConfidenceLikelihoodThreshold", doc = "for M2 events processed on each end, min threshold for CNQ of the M1 genotype of the end", optional = true)
    private Integer endConfidenceLikelihoodThreshold = null;

    @Argument(fullName = "processErrorProbs", doc = "process probabilities for CNV/quartet configurations that have an IBD error condition", optional = true)
    private Boolean processErrorProbs = false;

    @Argument(fullName = "m2Mode", doc = "mode for processing m2 events (doubleIBD or likelihood)", optional = true)
    private String m2Mode = SiblingIBDCopyNumberAnnotator.M2_MODE_DOUBLE_IBD;

    @Argument(fullName = "ends1VCF", doc = "VCF file with M1 genotypes of each end of M2 variants", optional = true)
    private File ends1VCF = null;

    @Argument(fullName = "endConfidenceReport", doc = "output file for report on end confidence for likelihood / M1 genotyping of M2 variants", optional = true)
    private File endConfidenceReport = null;

    @Argument(fullName="pedigree", shortName="ped", doc="Pedigree file")
    private File pedigreeFile = null;

    private SiblingIBDCopyNumberAnnotator annotator;

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();

        annotator = new SiblingIBDCopyNumberAnnotator();
        annotator.initialize(filterGenoypes,
                genotypeQualityThreshold,
                siblingIBDFile,
                minCallRate,
                endConfidenceLikelihoodThreshold,
                processErrorProbs,
                m2Mode,
                ends1VCF,
                endConfidenceReport,
                Collections.singletonList(pedigreeFile),
                true,
               out);
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        annotator.annotate(variant);
    }

    @Override
    public void closeTool() {
        super.closeTool();
        annotator.close();
    }
}
