package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class VcfToBedTest  extends TestUtils {

@Test(dataProvider="all-vcf-files")
public void test01(final String inputFile) 
	throws IOException
	{
	final File out = super.createTmpFile(".bed");
	Assert.assertEquals(0,new VcfToBed().instanceMain(new String[] {
		inputFile
		}));
	super.assertIsBed(out);
	}
@Test(dataProvider="all-vcf-files")
public void test02(final String inputFile) 
	throws IOException
	{
	final File out = super.createTmpFile(".bed");
	Assert.assertEquals(0,new VcfToBed().instanceMain(new String[] {
		"-header",
		inputFile
		}));
	super.assertIsBed(out);
	}

}
