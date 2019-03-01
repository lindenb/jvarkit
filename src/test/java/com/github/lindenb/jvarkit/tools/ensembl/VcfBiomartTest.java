package com.github.lindenb.jvarkit.tools.ensembl;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.stream.Collectors;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class VcfBiomartTest  {

	private final TestSupport support = new TestSupport();

	@DataProvider(name = "genomic-segments")
	public Object[][] createSomeGenomicSegments() {
		return new Object[][] {
			{"1",43_996_511,"CGCGAGCGCGAGGGGAGCGCGCGGCTGGAGCTGGCGCGGGAGCGGCGGGAGCGGTGGCGGCGGCAGAGGCGGCGGCTCCAGCTTCGGCTCC"},
			{"22",41_756_126,"CTAGAAAATAAACTTGTGCACTTTGACCTCTGTCCCCGAGATGT"},
			{"MT",5_881,"AGCCATTTTACCTCACCCCCACTGATGTTC"}
		};
		}
	private String xml_query = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" +
			"<!DOCTYPE Query>" +
			"<Query  virtualSchemaName = \"default\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.6\" >" +
			"\t\t\t" +
			"\t<Dataset name = \"hsapiens_gene_ensembl\" interface = \"default\" >" +
			"\t\t<Filter name = \"end\" value = \"10000000\"/>" +
			"\t\t<Filter name = \"start\" value = \"1000000\"/>" +
			"\t\t<Filter name = \"chromosome_name\" value = \"4\"/>" +
			"\t\t<Attribute name = \"ensembl_gene_id\" />" +
			"\t\t<Attribute name = \"ensembl_transcript_id\" />" +
			"\t\t<Attribute name = \"start_position\" />" +
			"\t\t<Attribute name = \"end_position\" />" +
			"\t\t<Attribute name = \"external_transcript_name\" />" +
			"\t\t<Attribute name = \"transcription_start_site\" />" +
			"\t\t<Attribute name = \"transcript_start\" />" +
			"\t\t<Attribute name = \"transcript_end\" />" +
			"\t</Dataset>" +
			"</Query>" 
			;

		
	@Test(dataProvider="genomic-segments",timeOut=180*1000L)
	public void test02(final String contig,int chromStart,String dna) 
		throws IOException
		{
		try {
			Path inVcfFile = support.createTmpPath(".vcf");
			PrintWriter pw = new PrintWriter(Files.newBufferedWriter(inVcfFile));
			pw.println("##fileformat=VCFv4.2");
			pw.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
			for(int i=0;i< dna.length();i++)
				{
				final String ref= String.valueOf(dna.charAt(i));
				pw.print(contig);
				pw.print("\t");
				pw.print(chromStart+i);
				pw.print("\t.\t");
				pw.print(ref);
				pw.print("\t");
				pw.print(Arrays.asList("A","C","G","T").stream().filter(A->!A.equals(ref)).collect(Collectors.joining(",")));
				pw.append("\t.\t.\t.");
				pw.println();
				}
			pw.flush();
			pw.close();
			support.assertIsVcf(inVcfFile);
			
			
			
			Path xmlQueryFile = support.createTmpPath(".xml");
			Writer xp = Files.newBufferedWriter(xmlQueryFile);
			xp.write(this.xml_query);
			xp.flush();
			xp.close();
			support.assertIsXml(xmlQueryFile);
			
			final Path out = support.createTmpPath(".vcf");
			final VcfBiomart cmd =new VcfBiomart();
			Assert.assertEquals(0,cmd.instanceMain(new String[] {
				"--contig","chromosome_name",
				"--start","start",
				"--end","end",
				"--xml",xmlQueryFile.toString(),
				"-o",out.toString(),
				inVcfFile.toString()
				}));
			Assert.assertEquals(support.variantStream(inVcfFile).count(),support.variantStream(out).count());
			} 
		finally
			{
			support.removeTmpFiles();
			}
		}

	
	
	
}
