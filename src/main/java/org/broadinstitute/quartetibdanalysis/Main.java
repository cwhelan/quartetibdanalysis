package org.broadinstitute.quartetibdanalysis;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Main extends org.broadinstitute.hellbender.Main {

    @Override
    protected String getCommandLineName() {
        return "QuartetIBDAnalysis";
    }

    @Override
    protected List<String> getPackageList() {
        final List<String> packageList = new ArrayList<>();
        packageList.addAll(Arrays.asList("org.broadinstitute.quartetibdanalysis"));
        return packageList;
    }

    public static void main(final String[] args) {
        new Main().mainEntry(args);
    }

}
