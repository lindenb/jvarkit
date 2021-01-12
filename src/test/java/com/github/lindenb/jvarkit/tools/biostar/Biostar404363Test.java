package com.github.lindenb.jvarkit.tools.biostar;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

import org.testng.Assert;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class Biostar404363Test {
	
private final TestSupport support = new TestSupport();
	
public void test01() throws IOException
	{
	try {
		final Path vcf = support.createTmpPath(".vcf");
		Files.write(vcf,Arrays.asList(
				"##fileformat=VCFv4.2" ,
				"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency among genotypes, for each ALT allele, in the same order as listed\">" ,
				"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" ,
				"14\t79838836\t.\tN\tA\t.\t.\tAF=0.5"
				));
		support.assertIsVcf(vcf);
		
		
		final Path out = support.createTmpPath(".bam");
		Assert.assertEquals(
				new Biostar404363().instanceMain(new String[] {
					"-o",out.toString(),
					"-vcf",vcf.toString(),
				support.resource("HG02260.transloc.chr9.14.bam")
				}),0);
		support.assertIsValidBam(out);
		}
	finally {
		support.removeTmpFiles();
		}
	}
}
