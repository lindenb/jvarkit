package com.github.lindenb.jvarkit.tools.vcfmerge;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.samtools.util.Interval;

public class VCFMergeTest extends TestUtils{


private File basetest(final String args) throws IOException {
	final File outvcf = super.createTmpFile(".vcf");
	final Object vcfpaths[]=new Object[5];
	for(int i=1;i<=vcfpaths.length;i++)
		{
		vcfpaths[i-1]=SRC_TEST_RESOURCE+"/S"+i+".vcf.gz";
		}
	Assert.assertEquals(new VCFMerge().instanceMain(newCmd().add(
			"-o",outvcf.getPath()).
			split(args).
			add(vcfpaths).
			make()),0);
	assertIsVcf(outvcf);
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
	final Interval interval = super.randomIntervalsFromDict(new File(SRC_TEST_RESOURCE+"/rotavirus_rf.fa"),1).get(0);
	basetest("--region "+interval.getContig()+":"+interval.getStart()+"-"+interval.getEnd());
	}
}
