package com.github.lindenb.jvarkit.tools.vcfmerge;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;


import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.samtools.util.Interval;


public class VCFMergeTest {


private Path basetest(final TestSupport support,final String args) throws IOException {
	
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
	final TestSupport support =new TestSupport();
	try {
		basetest(support,"");
		}
	finally {
		support.removeTmpFiles();
		}
	}


@Test
public void testRegion() throws IOException
	{
	final TestSupport support =new TestSupport();
	try {
	final Interval interval = support.randomIntervalsFromDict(Paths.get(support.resource("rotavirus_rf.fa")),1,1000).get(0);
		basetest(support,"--region "+interval.getContig()+":"+interval.getStart()+"-"+interval.getEnd());
		}
	finally {
		support.removeTmpFiles();
		}
	}
}
