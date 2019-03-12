package com.github.lindenb.jvarkit.tools.onekgenomes;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.samtools.util.IOUtil;

public class VcfAncestralAlleleTest {
	private final TestSupport support = new TestSupport();

	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][] {
			{support.resource("rotavirus_rf.fa"),support.resource("rotavirus_rf.vcf.gz")},
			{support.resource("toy.fa"),support.resource("toy.vcf.gz")}
		};
	}
	
	@Test(dataProvider="src1")
	public void test01(
			final String fasta,
			final String inputFile
			) throws IOException
	{
	try {
		final Path output = support.createTmpPath(".vcf");
		final Path manifest = support.createTmpPath(".MF");
		final File ref=new File(fasta);
		PrintWriter pw = new PrintWriter(Files.newBufferedWriter(manifest));
		for(final String l: IOUtil.slurpLines( new File(ref.getParentFile(),ref.getName()+".fai"))) {
			final String tokens[] = l.split("[\t]");
			pw.println(tokens[0]+"|"+tokens[0]+"xxx\t"+tokens[0]+"\t"+ref);
			}
		pw.flush();
		pw.close();
		 Assert.assertEquals(new VcfAncestralAllele().instanceMain(new String[] {
	     		"-o",output.toString(),
	     		"-m",manifest.toString(),
	     		inputFile
		 		}),0);
		support.assertIsVcf(output);
		
		Assert.assertTrue(support.variantStream(output).
				allMatch(V->V.hasAttribute("AA") && V.getReference().getDisplayString().equalsIgnoreCase(V.getAttributeAsString("AA", "."))));
		}
	finally
		{
		support.removeTmpFiles();
		}
	}
}
