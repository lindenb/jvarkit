package com.github.lindenb.jvarkit.tools.vcfstats;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.samtools.util.IOUtil;


public class VcfStatsTest  {
	

	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		final TestSupport support = new TestSupport();
		return support.toArrayArray(support.
				allVcfOrBcf().
				filter(F->F.endsWith(".vcf.gz")).
				map(F->new Object[] {F})
				)
				;
		}
	
	@Test(dataProvider="src1")
	public void basic(final String inputFile) 
		throws IOException
		{
		final TestSupport support = new TestSupport();
		try {
		final Path output = support.createTmpDirectory();
		List<String> args = new ArrayList<>();
		args.add("-o");
		args.add(output.toString());
		args.add(inputFile);
        Assert.assertEquals(0,new VcfStats().instanceMain(args));
        IOUtil.deleteDirectoryTree(output.toFile());
		} finally {
			support.removeTmpFiles();
		}
		}
	}
