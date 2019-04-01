package org.broadinstitute.quartetibdanalysis.tools;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;

public class QuartetIBDAnalysisProgramGroup implements CommandLineProgramGroup {

    @Override
    public String getName() { return "Quartet IBD Analysis Program Group"; }

    @Override
    public String getDescription() { return "Tools to analyze IBD state and copy number variation in quartets"; }

}
