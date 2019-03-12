package com.github.lindenb.jvarkit.tools.pcr;

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

import htsjdk.samtools.util.Interval;

@AlsoTest(LauncherTest.class)
public class PcrClipReadsTest {

	private final TestSupport support = new TestSupport();
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return support.toArrayArray(support.
				allSamOrBams().
				map(F->new Object[] {F})
				)
				;
		}
	
@Test(dataProvider="src1")
public void test01(final String bamPath) throws IOException {
	try {
	final Path out = support.createTmpPath(".bam");
	final Path bed = support.createTmpPath(".bed");
	
	for(final Interval rgn :support.randomIntervalsFromDict(Paths.get(bamPath),100,1_000)) {
		final PrintWriter pw = new PrintWriter(Files.newBufferedWriter(bed));
		pw.println(rgn.getContig()+"\t"+rgn.getStart()+"\t"+rgn.getEnd());
		pw.flush();
		pw.close();
		
		final String args[] = new String[] {
			"-o",out.toString(),
			"-B",bed.toString(),
			bamPath
		};
		
		final int ret = new PcrClipReads().instanceMain(args);
		if(ret!=0) {
			Assert.fail(String.join(" ",args)+" "+rgn+" "+ret);
			}
		
		support.assertIsValidBam(out);
		}
	} finally {
		support.removeTmpFiles();
	}
	}
	
}
