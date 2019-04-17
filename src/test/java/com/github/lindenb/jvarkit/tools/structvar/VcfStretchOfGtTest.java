package com.github.lindenb.jvarkit.tools.structvar;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.Reporter;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class VcfStretchOfGtTest{
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
	public void test01(String inputFile) throws IOException
		{
		try {
		final Path ped = support.createRandomPedigreeFromFile(inputFile);
		if(ped==null) {
			 Reporter.getCurrentTestResult().setAttribute("warn", "No Pedigree for "+inputFile);
			return;
			}

		final Path out = support.createTmpPath(".tsv");
		Assert.assertEquals(new VcfStretchOfGt().instanceMain(new String[] {
				"-o",out.toString(),
				"--pedigree",
				ped.toString(),
				inputFile
				}),0
			);
		support.assertTsvTableIsConsitent(out, null);
		} finally {
			support.removeTmpFiles();
		}
		}
	}
