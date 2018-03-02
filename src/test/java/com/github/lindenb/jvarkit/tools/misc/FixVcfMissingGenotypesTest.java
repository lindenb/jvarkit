package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.util.stream.Collectors;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;


public class FixVcfMissingGenotypesTest extends TestUtils
	{
	@Test
	public void test01() throws Exception{
		File nocallvcffile = super.createTmpFile(".vcf");
		VCFFileReader r1=new VCFFileReader(
				new File("./src/test/resources/rotavirus_rf.vcf.gz"),
				false);
		final VariantContextWriter vcw = VCFUtils.createVariantContextWriter(nocallvcffile);
		vcw.writeHeader(r1.getFileHeader());
		r1.iterator().stream().forEach(V->{
			VariantContextBuilder vcb = new VariantContextBuilder(V);
			
			vcb.genotypes(V.getGenotypes().stream().map(G->GenotypeBuilder.
						createMissing(G.getSampleName(), 2)).
						collect(Collectors.toList()));
			vcw.add(vcb.make());
			});
		vcw.close();
		r1.close();
		final File output = super.createTmpFile(".vcf");

		final FixVcfMissingGenotypes cmd =new FixVcfMissingGenotypes();
		Assert.assertEquals(0,cmd.instanceMain(newCmd().add(
			"-B","./src/test/resources/S1.bam", 
			"-B","./src/test/resources/S2.bam", 
			"-B","./src/test/resources/S3.bam", 
			"-B","./src/test/resources/S4.bam", 
			"-B","./src/test/resources/S5.bam", 
			"-o",output.getPath(),
			nocallvcffile
			).make()
			));
		assertIsVcf(output);
		Assert.assertFalse(
				super.variantStream(output).
				flatMap(V->V.getGenotypes().stream()).
				noneMatch(G->G.isHomRef()),
				"at least one GT should be 0/0"
				);
		}
	}