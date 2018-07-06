package com.github.lindenb.jvarkit.tools.xcontamination;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class XContaminationsTest extends TestUtils {
public void test01() throws IOException {
	
	File output = super.createTmpFile(".vcf");
	Assert.assertEquals(new XContaminations().instanceMain(
    		newCmd().add(
    		"-ov","-sample",
    		"-o",output,
    		SRC_TEST_RESOURCE+"/rotavirus_rf.vcf.gz",
    		SRC_TEST_RESOURCE+"/S1.bam",
    		SRC_TEST_RESOURCE+"/S2.bam",
    		SRC_TEST_RESOURCE+"/S3.bam",
    		SRC_TEST_RESOURCE+"/S4.bam",
    		SRC_TEST_RESOURCE+"/S5.bam"
    			).make()
    	),0);
	super.assertIsVcf(output);
	}
}
