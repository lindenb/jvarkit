package com.github.lindenb.jvarkit.tools.structvar.gridss;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;
import com.github.lindenb.jvarkit.variant.vcf.BcfToolsReaderTest;

@AlsoTest({LauncherTest.class,BcfToolsReaderTest.class})
public class GridssMergeBndTest {
	private final TestSupport support = new TestSupport();
	@Test
	public void test01() throws IOException {
		try {
			final String header=
				"##fileformat=VCFv4.2\n"+
				"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"+
				"##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"CI\">\n"+
				"##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n"+
				"##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">\n"+
				"##contig=<ID=chr1,length=249250621,assembly=HS37D5>\n"
				;
			Path vcf1 = support.createTmpPath(".vcf");
			PrintWriter pw = new PrintWriter(Files.newBufferedWriter(vcf1));
			pw.print(header);
			pw.print(
				"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" +
				"chr1\t10285\t.\tT\t]chr1:10248]ACCCTAACCCTAACCCTAACCCCTAACCCTAACCCCT\t156.28\t.\t.\tGT\t.\n" +
				"chr1\t54783\t.\tT\t]chr1:54777]CTTTCT\t223.47\t.\t.\tGT\t.\n" 	
				);
			pw.flush();
			pw.close();
			support.assertIsVcf(vcf1);
			Path vcf2 = support.createTmpPath(".vcf");
			pw = new PrintWriter(Files.newBufferedWriter(vcf2));
			pw.print(header);
			pw.print(
				"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2\n" +
				"chr1\t54783\t.\tT\t]chr1:54777]CTTTCT\t223.47\t.\t.\tGT\t.\n" 	+
				"chr1\t55285\t.\tT\t]chr1:10248]ACCCTAACCCTAACCCTAACCCCTAACCCTAACCCCT\t156.28\t.\t.\tGT\t.\n" 
				);
			pw.flush();
			pw.close();
			support.assertIsVcf(vcf2);
			Path vcf3 = support.createTmpPath(".vcf");

			
			Assert.assertEquals(
					new GridssMergeBnd().instanceMain(new String[]{
					"-o",vcf3.toString(),
					vcf1.toString(),
					vcf2.toString()
					})
					, 0);
			support.assertIsVcf(vcf3);
			}	
		finally {
			support.removeTmpFiles();
		}
	}
}
