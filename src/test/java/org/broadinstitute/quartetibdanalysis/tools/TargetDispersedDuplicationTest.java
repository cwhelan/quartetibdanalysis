package org.broadinstitute.quartetibdanalysis.tools;

import org.testng.Assert;
import org.testng.annotations.Test;

public class TargetDispersedDuplicationTest {

    @Test
    public void testCalcPwrong() throws Exception {
        Assert.assertEquals(new Double(-1 * Math.log10(TargetDispersedDuplication.QuartetCopyNumberEvent.calcPWrong(99, 90.5, 81.2, 99))).intValue() , 8);
    }

}