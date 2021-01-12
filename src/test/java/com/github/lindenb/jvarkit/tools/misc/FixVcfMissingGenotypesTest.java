package com.github.lindenb.jvarkit.tools.misc;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.stream.Collectors;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VCFUtilsTest;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFReader;

@AlsoTest(VCFUtilsTest.class)
public class FixVcfMissingGenotypesTest 
	{
	@Test
	public void test01() throws Exception{
		 final TestSupport support = new TestSupport();
			 try {
			Path nocallvcffile = support.createTmpPath(".vcf");
			VCFReader r1= VCFReaderFactory.makeDefault().open(
					Paths.get(support.resource("rotavirus_rf.vcf.gz")),
					false);
			final VariantContextWriter vcw = VCFUtils.createVariantContextWriter(nocallvcffile.toFile());
			vcw.writeHeader(r1.getHeader());
			r1.iterator().stream().forEach(V->{
				VariantContextBuilder vcb = new VariantContextBuilder(V);
				
				vcb.genotypes(V.getGenotypes().stream().map(G->GenotypeBuilder.
							createMissing(G.getSampleName(), 2)).
							collect(Collectors.toList()));
				vcw.add(vcb.make());
				});
			vcw.close();
			r1.close();
			final Path output = support.createTmpPath(".vcf");
	
			final FixVcfMissingGenotypes cmd =new FixVcfMissingGenotypes();
			Assert.assertEquals(0,cmd.instanceMain(new String[] {
				"-B",support.resource("S1.bam"), 
				"-B",support.resource("S2.bam"), 
				"-B",support.resource("S3.bam"), 
				"-B",support.resource("S4.bam"), 
				"-B",support.resource("S5.bam"), 
				"-o",output.toString(),
				nocallvcffile.toString()
				}));
			support.assertIsVcf(output);
			Assert.assertFalse(
					support.variantStream(output).
					flatMap(V->V.getGenotypes().stream()).
					noneMatch(G->G.isHomRef()),
					"at least one GT should be 0/0"
					);
			 }
		finally {
			 support.removeTmpFiles();
		 }
		}
	}