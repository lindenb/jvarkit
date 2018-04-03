package com.github.lindenb.jvarkit.tools.pcr;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.samtools.util.Interval;

public class PcrClipReadsTest extends TestUtils {

@Test(dataProvider="all-sam-or-bam-files")
public void test01(final String bamPath) throws IOException {
	final File out = super.createTmpFile(".bam");
	final File bed = super.createTmpFile(".bed");
	
	for(final Interval rgn :randomIntervalsFromDict(new File(bamPath),100)) {
		final PrintWriter pw = new PrintWriter(bed);
		pw.println(rgn.getContig()+"\t"+rgn.getStart()+"\t"+rgn.getEnd());
		pw.flush();
		pw.close();
		
		final String args[] = newCmd().
			add("-o",out.getPath()).
			add("-B",bed.getPath()).
			add(bamPath).make()
			;
		
		final int ret = new PcrClipReads().instanceMain(args);
		if(ret!=0) {
			Assert.fail(String.join(" ",args)+" "+rgn+" "+ret);
			}
		
		assertIsValidBam(out);
		}
	}
	
}
