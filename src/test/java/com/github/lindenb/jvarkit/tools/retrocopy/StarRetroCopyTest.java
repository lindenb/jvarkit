package com.github.lindenb.jvarkit.tools.retrocopy;


import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

import org.testng.Assert;
import org.testng.annotations.Test;

@AlsoTest(LauncherTest.class)
public class StarRetroCopyTest   {
	private final TestSupport support = new TestSupport();
	@Test
	public void test01() throws IOException {
		try {
			final Path gtf = support.createTmpPath(".gtf");
			final Path bedpe = support.createTmpPath(".bedpe");
			Files.write(gtf,Arrays.asList(
					"4\tensembl\ttranscript\t15606631\t15656964\t.\t-\t.\tgene_id \"ENSG00000118564\"; gene_version \"10\"; transcript_id \"ENST00000382358\"; transcript_version \"4\"; gene_name \"FBXL5\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"FBXL5-201\"; transcript_source \"ensembl\"; transcript_biotype \"protein_coding\"; tag \"basic\";" +
					"4\tensembl\texon\t15656827\t15656964\t.\t-\t.\tgene_id \"ENSG00000118564\"; gene_version \"10\"; transcript_id \"ENST00000382358\"; transcript_version \"4\"; exon_number \"1\"; gene_name \"FBXL5\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"FBXL5-201\"; transcript_source \"ensembl\"; transcript_biotype \"protein_coding\"; exon_id \"ENSE00003576135\"; exon_version \"1\"; tag \"basic\";" +
					"4\tensembl\texon\t15646116\t15646331\t.\t-\t.\tgene_id \"ENSG00000118564\"; gene_version \"10\"; transcript_id \"ENST00000382358\"; transcript_version \"4\"; exon_number \"2\"; gene_name \"FBXL5\"; gene_source \"ensembl_havana\"; gene_biotype \"protein_coding\"; transcript_name \"FBXL5-201\"; transcript_source \"ensembl\"; transcript_biotype \"protein_coding\"; exon_id \"ENSE00003575511\"; exon_version \"1\"; tag \"basic\";" 
					));
			
			final Path out = support.createTmpPath(".vcf");
			final Path out2 = support.createTmpPath(".bam");
			Assert.assertEquals(new StarRetroCopy().instanceMain(new String[] {
					"--gtf",gtf.toString(),
					"-o",out.toString(),
					"--bam",out2.toString(),
					"--bedpe",bedpe.toString(),
					"--bwa",
					support.resource("retrocopy01.bwa.bam")
					}),0);
			
			support.assertIsVcf(out);
			support.assertIsValidBam(out2);
		} finally
		{
			support.removeTmpFiles();
		}
	}
}
