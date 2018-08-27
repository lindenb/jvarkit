package com.github.lindenb.jvarkit.tools.sam4weblogo;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.samtools.util.Interval;


public class SAM4WebLogoTest extends TestUtils
	{
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new ParamCombiner().
			initList(collectIndexedBams()).
			build();
		}
	private File basetest(final String inBam,String params,boolean fastq) throws IOException {
		final File out = createTmpFile(fastq?".fastq":".txt");
		for(final Interval interval: super.randomIntervalsFromDict(new File(inBam), 20)) {
			if(interval.getContig().contains(":")) continue;
			final String args[]=newCmd().add(
					"-o",out.getPath(),
					"-r",interval.getContig()+":"+interval.getStart()+"-"+interval.getEnd()).
					split(params).
					add(inBam).
					make();
			
			final int ret=new SAM4WebLogo().instanceMain(args);
			if(ret!=0) {
				Assert.fail(Arrays.asList(args).toString());
				}
			}
		return out;
		}
	
	@Test(dataProvider="src1")
	public void testNoClip(final String inBam) throws IOException {
		basetest(inBam,"",false);
		}
	@Test(dataProvider="src1")
	public void testClip(final String inBam) throws IOException {
		basetest(inBam,"-c",false);
		}
	@Test(dataProvider="src1")
	public void testFastq(final String inBam) throws IOException {
		final File out=basetest(inBam,"-c --fastq -fqu _ -fqp _",true);
		assertIsFastq(out);
		}
	}
