package com.github.lindenb.jvarkit.tools.vcfmerge;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

import htsjdk.samtools.util.Interval;

@AlsoTest(LauncherTest.class)
public class VCFMergeTest {

	private final TestSupport support =new TestSupport();

private Path basetest(final String args) throws IOException {
	final Path outvcf = support.createTmpPath(".vcf");
	
	List<String> al = new ArrayList<>();
	al.add("-o");
	al.add(outvcf.toString());
	
	Arrays.stream(args.split("[ \t]")).filter(S->!S.isEmpty()).forEach(S->al.add(S));
	
	for(int i=1;i<=5;i++)
		{
		al.add(support.resource("S"+i+".vcf.gz"));
		}
	
	
	Assert.assertEquals(new VCFMerge().instanceMain(al),0);
	support.assertIsVcf(outvcf);
	return outvcf;
	}

@Test
public void testUnsorted() throws IOException
	{
	basetest("");
	}

@Test
public void testSorted() throws IOException
	{
	basetest("--sorted");
	}

@Test
public void testRegion() throws IOException
	{
	final Interval interval = support.randomIntervalsFromDict(Paths.get(support.resource("rotavirus_rf.fa")),1,1000).get(0);
	basetest("--region "+interval.getContig()+":"+interval.getStart()+"-"+interval.getEnd());
	}
}
