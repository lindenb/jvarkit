package com.github.lindenb.jvarkit.tools.vcfstats;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class VcfStatsTest  {
	
	private final TestSupport support = new TestSupport();

	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return support.toArrayArray(support.
				allVcfOrBcf().
				map(F->new Object[] {F})
				)
				;
		}
	
	@Test(dataProvider="src1")
	public void test01(final String inputFile) 
		throws IOException
		{
		try {
		final Path output = support.createTmpPath(".zip");
		final Path ped = support.createRandomPedigreeFromFile(inputFile);
		List<String> args = new ArrayList<>();
		args.add("-o");
		args.add(output.toString());
		if(ped!=null) {
			args.add("--pedigree");
			args.add(ped.toString());
		}
		args.add(inputFile);
		
        Assert.assertEquals(0,new VcfStats().instanceMain(args));
        support.assertZip(output);
		} finally {
			support.removeTmpFiles();
		}
		}
	}
