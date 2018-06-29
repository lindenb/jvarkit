package com.github.lindenb.jvarkit.tools.vcfcomposite;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class VCFCompositeTest extends TestUtils{
	@Test(dataProvider="all-vcf-files")
	public void test01(final String inputFile) 
		throws IOException
		{
		final File ped = super.createRandomPedigreeFromFile(inputFile);
		if(ped==null) return;
		final File output = super.createTmpFile(".tsv");
        Assert.assertEquals(new VCFComposite().instanceMain(
        		newCmd().add(
        		"-o",output.getPath(),
        		"-m","RecessiveComposite",
        		"--pedigree",ped).
        		add(inputFile).make()
        	),0);
        super.assertTsvTableIsConsitent(output, null);
		}
}
