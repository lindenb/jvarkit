package com.github.lindenb.jvarkit.tools.bioalcidae;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;
import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class BioAlcidaeTest extends TestUtils
	{
	
	@Test(dataProvider = "all-vcf-files")
	public void testVCF(final String vcfpath) 
		throws IOException
		{
		final File out = super.createTmpFile(".txt");
		
		Assert.assertEquals(new BioAlcidae().instanceMain(newCmd().
				add("-o",out.getPath()).
				add("-e","while(iter.hasNext()) out.println(iter.next().getContig());").
				add(vcfpath).make()
				),0);
		assertIsNotEmpty(out);
		}
	
}
