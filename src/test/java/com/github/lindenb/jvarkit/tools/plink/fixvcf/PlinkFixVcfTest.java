package com.github.lindenb.jvarkit.tools.plink.fixvcf;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.bio.AcidNucleics;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class PlinkFixVcfTest {

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][]{
			{"rotavirus_rf.vcf.gz","rotavirus_rf.fa"},
			};
		}
	@Test(dataProvider="src1")
	public void testFixVcf(final String vcf,final String ref) throws IOException {
		final TestSupport support = new TestSupport();
		try {
			int n_variants1=0;
			final Path tmp = support.createTmpPath(".vcf");
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
							String[] tokens = CharSplitter.TAB.split(line);
							if(!AcidNucleics.isATGC(tokens[4])) continue;//indel
							String swap = tokens[3];
							tokens[3]=tokens[4];
							tokens[4]=swap;
							pw.println(String.join("\t", tokens));
							n_variants1++;
							}
						}
					}
				pw.flush();
				}
			
			
			final Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(
				new PlinkFixVcf().instanceMain(new String[] {
				"-R",support.resource(ref),
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
