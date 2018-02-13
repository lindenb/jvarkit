package com.github.lindenb.jvarkit.tools.ensembl;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.stream.Collectors;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class VcfEnsemblVepRestTest extends TestUtils {
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][] {
			{"1 44044896 C G"}
		};
	}
	private String toVcfLine(final String vcfline) {
		String tokens[]= vcfline.split("[ \t]+");
		Assert.assertEquals(tokens.length, 4);
		StringBuilder out=new StringBuilder();
		out.append(tokens[0]);
		out.append("\t");
		out.append(tokens[1]);
		out.append("\t");
		out.append(".");
		out.append("\t");
		out.append(tokens[2]);
		out.append("\t");
		out.append(tokens[3]);
		out.append("\t.\t.\t.");
		return out.toString();
	}
	
	@Test(dataProvider="src1")
	public void test01(final String vcfline) 
		throws IOException
		{
		File inVcfFile = super.createTmpFile(".vcf");
		PrintWriter pw = new PrintWriter(inVcfFile);
		pw.println("##fileformat=VCFv4.2");
		pw.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
		pw.println(toVcfLine(vcfline));
		pw.flush();
		pw.close();
		
		
		final File out = super.createTmpFile(".vcf");
		final VcfEnsemblVepRest cmd =new VcfEnsemblVepRest();
		Assert.assertEquals(0,cmd.instanceMain(new String[] {
			"-o",out.getPath(),
			inVcfFile.getPath()
			}));
		
		}
	
	@DataProvider(name = "src2")
	public Object[][] createData2() {
		return new Object[][] {
			{"1",43_996_511,"CGCGAGCGCGAGGGGAGCGCGCGGCTGGAGCTGGCGCGGGAGCGGCGGGAGCGGTGGCGGCGGCAGAGGCGGCGGCTCCAGCTTCGGCTCC"},
			{"22",41_756_126,"CTAGAAAATAAACTTGTGCACTTTGACCTCTGTCCCCGAGATGT"},
			{"MT",5_881,"AGCCATTTTACCTCACCCCCACTGATGTTC"}
		};
		}
		
	@Test(dataProvider="src2")
	public void test02(final String contig,int chromStart,String dna) 
		throws IOException
		{
		File inVcfFile = super.createTmpFile(".vcf");
		PrintWriter pw = new PrintWriter(inVcfFile);
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
			pw.print(ref+"A");//insertion
			pw.append("\t.\t.\t.");
			pw.println();
			}
		pw.flush();
		pw.close();
		
		
		final File out = super.createTmpFile(".vcf");
		final VcfEnsemblVepRest cmd =new VcfEnsemblVepRest();
		Assert.assertEquals(0,cmd.instanceMain(new String[] {
			"-o",out.getPath(),
			inVcfFile.getPath()
			}));
		
		}

	
	
	
}
