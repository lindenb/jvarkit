package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;


public class VcfRemoveUnusedAltTest extends TestUtils {
	@Test(dataProvider="all-vcf-files")
	public void test01(final String inputFile) 
		throws IOException
		{		
		
		final File output = super.createTmpFile(".vcf");
        Assert.assertEquals(0,new VcfRemoveUnusedAlt().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"--onespan",
        		inputFile
        	}));
        assertIsVcf(output);
		}
	}
