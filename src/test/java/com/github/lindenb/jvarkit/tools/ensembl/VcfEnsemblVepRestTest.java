package com.github.lindenb.jvarkit.tools.ensembl;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.stream.Collectors;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class VcfEnsemblVepRestTest {
	
	private final TestSupport support = new TestSupport();

	
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
		try {
		Path inVcfFile = support.createTmpPath(".vcf");
		PrintWriter pw = new PrintWriter(Files.newBufferedWriter(inVcfFile));
		pw.println("##fileformat=VCFv4.2");
		pw.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
		pw.println(toVcfLine(vcfline));
		pw.flush();
		pw.close();
		
		
		final Path out = support.createTmpPath(".vcf");
		final VcfEnsemblVepRest cmd =new VcfEnsemblVepRest();
		Assert.assertEquals(0,cmd.instanceMain(new String[] {
			"-o",out.toString(),
			inVcfFile.toString()
			}));
		} finally {
			support.removeTmpFiles();
		}
		}
	
		
	@Test(dataProvider="genomic-segments",timeOut=180*1000)
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
			pw.print(ref+"A");//insertion
			pw.append("\t.\t.\t.");
			pw.println();
			}
		pw.flush();
		pw.close();
		
		
		final Path out = support.createTmpPath(".vcf");
		final VcfEnsemblVepRest cmd =new VcfEnsemblVepRest();
		Assert.assertEquals(0,cmd.instanceMain(new String[] {
			"-o",out.toString(),
			inVcfFile.toString()
			}));
		Assert.assertEquals(support.variantStream(inVcfFile).count(),support.variantStream(out).count());
		} finally {
			support.removeTmpFiles();
		}
		}

	
	
	
}
