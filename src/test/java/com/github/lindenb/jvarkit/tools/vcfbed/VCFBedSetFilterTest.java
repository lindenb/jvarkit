package com.github.lindenb.jvarkit.tools.vcfbed;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.testng.Assert;
import org.testng.annotations.Test;
import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class VCFBedSetFilterTest extends TestUtils{
	
	@Test(dataProvider = "all-vcf-files")
	public void testMemory(final String invcf) throws IOException {
		final File out = createTmpFile(".vcf");
		final File bedout = createTmpFile(".bed");
		
		final PrintWriter pw = new PrintWriter(bedout);
		super.randomIntervalsFromDict(new File(invcf), 100).stream().forEach(
				R->pw.println(R.getContig()+"\t"+(R.getStart()-1)+"\t"+R.getEnd())
				);
		pw.flush();
		pw.close();
		
		
		Assert.assertEquals(new VCFBedSetFilter().instanceMain(new String[] {
			"-o",out.getPath(),
			"-m",bedout.getPath(),
			invcf
			}),0);
		assertIsVcf(out);
		}

}
