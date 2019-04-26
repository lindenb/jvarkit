package com.github.lindenb.jvarkit.tools.vcfbed;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class VCFBedSetFilterTest {
	private final TestSupport support =new TestSupport();

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return support.toArrayArray(support.
				allVcfOrBcf().
				map(F->new Object[] {F})
				)
				;
		}
	
	@Test(dataProvider = "src1")
	public void testMemory(final String invcf) throws IOException {
		try {
		final Path out = support.createTmpPath(".vcf");
		final Path bedout = support.createTmpPath(".bed");
		
		final PrintWriter pw = new PrintWriter(Files.newBufferedWriter(bedout));
		support.randomIntervalsFromDict(Paths.get(invcf), 100,1000).stream().forEach(
				R->pw.println(R.getContig()+"\t"+(R.getStart()-1)+"\t"+R.getEnd())
				);
		pw.flush();
		pw.close();
		
		
		Assert.assertEquals(new VCFBedSetFilter().instanceMain(new String[] {
			"-o",out.toString(),
			"--memory",
			"-m",bedout.toString(),
			invcf
			}),0);
		support.assertIsVcf(out);
		} finally {
			support.removeTmpFiles();
			}
		}

}
