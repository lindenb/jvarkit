package com.github.lindenb.jvarkit.tools.vcfsetdict;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.dict.SequenceDictionaryExtractor;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.tools.vcfrebase.VcfRebase;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class VcfSetSequenceDictionaryTest {

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][]{
			{"rotavirus_rf.vcf.gz","rotavirus_rf.fa"},
			};
		}
	@Test(dataProvider="src1")
	public void test1(final String vcf,final String ref) throws IOException {
		final TestSupport support = new TestSupport();
		try {
			int n_variants1=0;
			Path tmp = support.createTmpPath(".vcf");
			try(PrintWriter pw = IOUtils.openPathForPrintWriter(tmp)) {
			try(BufferedReader r=IOUtils.openPathForBufferedReading(Paths.get(support.resource(vcf)))) {
					for(;;) {
						String line = r.readLine();
						if(line==null) break;
						if(line.startsWith("##contig=")) {
							continue;
							}
						else if(line.startsWith("#")) {
							pw.println(line);
							}
						else
							{
							pw.println("chr"+line);
							n_variants1++;
							}
						}
					}
				pw.flush();
				}
			
			final Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(
				new VcfRebase().instanceMain(new String[] {
				"-R",support.resource(ref),
				"-n","SKIP",
				"-o",out.toString(),
				tmp.toString()
				}),0);
			support.assertIsVcf(out);
			final SAMSequenceDictionary dict= SAMSequenceDictionaryExtractor.extractDictionary(Paths.get(support.resource(ref)));
			SequenceUtil.assertSequenceDictionariesEqual(
				SAMSequenceDictionaryExtractor.extractDictionary(out),
				dict
				);
			support.variantStream(out)
				.map(V->V.getContig())
				.forEach(C->Assert.assertNotNull(dict.getSequence(C),String.valueOf(C)) );
			Assert.assertEquals(support.variantStream(out).count(), n_variants1);
			}
		finally {
			support.removeTmpFiles();
			}
		}

}
