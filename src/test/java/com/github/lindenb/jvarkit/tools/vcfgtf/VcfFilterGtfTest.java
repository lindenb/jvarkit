package com.github.lindenb.jvarkit.tools.vcfgtf;

import java.nio.file.Path;
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
			final Path outVcf = support.createTmpPath(".vcf");
			final List<String> args= new ArrayList<>();
			args.add("-o");
			args.add(outVcf.toString());
			args.add("--gtf");
			args.add(support.resource("Homo_sapiens.GRCh37.87.gtf.gz"));
			args.add(support.resource("gnomad.exomes.r2.0.1.sites.vcf.gz"));
			
			Assert.assertEquals(new VcfFilterGtf().instanceMain(args), 0);
		}
		finally {
			support.removeTmpFiles();
		}
		
	}
}
