package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.io.PrintWriter;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.samtools.util.Interval;

import org.testng.Assert;
import org.testng.annotations.Test;

public class BedNonOverlappingSetTest extends TestUtils {
	@Test
    public void test01() throws IOException {
	final File tmpBed =super.createTmpFile(".bed");
    final PrintWriter pw = new PrintWriter(tmpBed);
	for(final Interval i:super.randomIntervalsFromDict(new File(SRC_TEST_RESOURCE+"/rotavirus_rf.fa"), 1000))
		{
		pw.println(i.getContig()+"\t"+(i.getStart()-1)+"\t"+i.getEnd());
		}
	pw.flush();
	pw.close();
	assertIsBed(tmpBed);
	File pat= File.createTempFile("tmp.__SETID__.", ".bed");
	Assert.assertEquals(
			new BedNonOverlappingSet().instanceMain(newCmd().add(
					"-o",pat.getPath(),
					"-x","10",
					"-R",SRC_TEST_RESOURCE+"/rotavirus_rf.fa",
					tmpBed.getPath()
			).make()),0);
	for(final File bf:pat.getParentFile().listFiles(new FileFilter() {
		@Override
		public boolean accept(File f) {
			return f.isFile() && f.exists() && f.getName().startsWith("tmp.") && f.getName().endsWith(".bed");
			}
		}))
		{
		assertIsBed(bf);
		bf.delete();
		}
	
	pat.delete();
	}
}
