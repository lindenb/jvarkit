package com.github.lindenb.jvarkit.tools.vcfgtf;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class VcfFilterGtfTest {
	@Test
	public void doBasicTest() throws Exception {
		final  TestSupport support =new TestSupport();
		try {
			final Path bed = support.createTmpBed(Paths.get(support.resource("rotavirus_rf.dict")));
			final Path outBed = support.createTmpPath(".bed");
			final List<String> args= new ArrayList<>();
			args.add("-o");
			args.add(outBed.toString());
			args.add("-R");
			args.add(support.resource("rotavirus_rf.dict"));
			args.add("--gtf");
			args.add(support.resource("rotavirus_rf.gtf.gz"));
			args.add(bed.toString());
			
			Assert.assertEquals(new VcfFilterGtf().instanceMain(args), 0);
			support.assertIsBed(outBed);
		}
		finally {
			support.removeTmpFiles();
		}
		
	}
}
