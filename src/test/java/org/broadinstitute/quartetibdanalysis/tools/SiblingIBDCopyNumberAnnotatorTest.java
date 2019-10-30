package org.broadinstitute.quartetibdanalysis.tools;

import org.testng.Assert;
import org.testng.annotations.Test;

public class SiblingIBDCopyNumberAnnotatorTest {

    @Test
    public void testIBD0Consistent() throws Exception {
        final SiblingIBDCopyNumberAnnotator annotator = new SiblingIBDCopyNumberAnnotator();
        Assert.assertTrue(annotator.consistentWithIDB0(3, 2, 3, 2));
        Assert.assertFalse(annotator.consistentWithIDB0(3, 2, 2, 2));
    }

    @Test
    public void testIBD1PConsistent() throws Exception {
        final SiblingIBDCopyNumberAnnotator annotator = new SiblingIBDCopyNumberAnnotator();
        Assert.assertTrue(annotator.consistentWithIDB1P(2, 2, 3, 2));
        Assert.assertFalse(annotator.consistentWithIDB1P(1, 1, 2, 1));
    }

    @Test
    public void testIBD2Consistent() throws Exception {
        final SiblingIBDCopyNumberAnnotator annotator = new SiblingIBDCopyNumberAnnotator();
        Assert.assertTrue(annotator.consistentWithIDB2(3, 3));
        Assert.assertFalse(annotator.consistentWithIDB2(3, 2));
    }

}