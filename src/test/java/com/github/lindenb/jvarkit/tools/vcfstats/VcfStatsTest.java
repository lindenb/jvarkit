package com.github.lindenb.jvarkit.tools.vcfstats;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;


public class VcfStatsTest extends TestUtils {
	@Test(dataProvider="all-vcf-files")
	public void test01(final String inputFile) 
		throws IOException
		{
		final File output = super.createTmpFile(".zip");
		final File ped = super.createRandomPedigreeFromFile(inputFile);
        Assert.assertEquals(0,new VcfStats().instanceMain(
        		newCmd().add(
        		"-o",output.getPath()).
        		addIf(ped!=null, "--pedigree",ped).
        		add(inputFile).make()
        	));
		}
	}
