package com.github.lindenb.jvarkit.tools.vcfcomposite;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

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
		if(Files.lines(ped.toPath()).noneMatch(L->L.endsWith("1")))
			{
			//no affected in pedigree ?
			return;
			}
		if(Files.lines(ped.toPath()).noneMatch(L->L.endsWith("0")))
			{
			//no unaffected in pedigree ?
			return;
			}
		final File output = super.createTmpFile(".vcf");
        Assert.assertEquals(new VCFComposite().instanceMain(
        		newCmd().add(
        		"-o",output.getPath(),
        		"--pedigree",ped).
        		add(inputFile).make()
        	),0);
        super.assertIsVcf(output);
		}
}
