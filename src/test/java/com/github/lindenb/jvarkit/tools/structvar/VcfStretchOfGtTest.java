package com.github.lindenb.jvarkit.tools.structvar;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.testng.Assert;
import org.testng.Reporter;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class VcfStretchOfGtTest  extends TestUtils{
	@Test(dataProvider="all-vcf-files")
	public void test01(String inputFile) throws IOException
		{
		final File ped = super.createRandomPedigreeFromFile(inputFile);
		if(ped==null) {
			 Reporter.getCurrentTestResult().setAttribute("warn", "No Pedigree for "+inputFile);
			return;
			}

		final File out = super.createTmpFile(".tsv");
		Assert.assertEquals(new VcfStretchOfGt().instanceMain(new String[] {
				"-o",out.getPath(),
				"--pedigree",
				ped.getPath(),
				inputFile
				}),0
			);
		assertTsvTableIsConsitent(out, null);
		}
	}
