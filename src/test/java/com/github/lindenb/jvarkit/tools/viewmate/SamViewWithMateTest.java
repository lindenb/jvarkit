package com.github.lindenb.jvarkit.tools.viewmate;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class SamViewWithMateTest  extends TestUtils {
@Test
public void test01() throws IOException {
	final File input = new File(SRC_TEST_RESOURCE+"/HG02260.transloc.chr9.14.bam");
	final File output = super.createTmpFile(".bam");
	Assert.assertEquals(0,new SamViewWithMate().instanceMain(
    		newCmd().add(
    		"-o",output,
    		"-r","9:137230721-137230796",
    		input
    		).make()
    	));
    assertIsValidBam(output);
    SamReaderFactory srf=SamReaderFactory.makeDefault();
	for(final String str: new String[]{"9","14"}) {
		try(SamReader sr = srf.open(output)) {
			Assert.assertTrue(sr.iterator().stream().anyMatch(R->str.equals(R.getContig())));
			}
		}
    
	}
}
