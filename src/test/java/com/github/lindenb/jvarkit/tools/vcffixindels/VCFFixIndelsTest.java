package com.github.lindenb.jvarkit.tools.vcffixindels;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class VCFFixIndelsTest {
	
	private final TestSupport support =new TestSupport();

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
		final Path output = support.createTmpPath(".vcf");
		Assert.assertEquals(new VCFFixIndels().instanceMain(new String[] {
        		"-o",output.toString(),
        		inputFile.toString()
        		}
        	),0);
		support.assertIsVcf(output);
		}
	}
