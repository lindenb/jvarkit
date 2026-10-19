package com.github.lindenb.jvarkit.tools.vcfsetdict;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;

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
	public void testNoDict(final String vcf,final String ref) throws IOException {
		run_test(0,false,vcf,ref);
		}
	@Test(dataProvider="src1")
	public void testDictChr(final String vcf,final String ref) throws IOException {
		run_test(1,false,vcf,ref);
		}
	
	@Test(dataProvider="src1")
	public void testPlainNoDict(final String vcf,final String ref) throws IOException {
		run_test(0,true,vcf,ref);
		}
	@Test(dataProvider="src1")
	public void testPlainDict(final String vcf,final String ref) throws IOException {
		run_test(1,true,vcf,ref);
		}

	private void run_test(int type,boolean plain,final String vcf,final String ref) throws IOException {

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
							if(type==0) {
								continue;
								}
							else if(type==1) {
								pw.println(line.replaceAll("ID=", "ID=chr"));
								}
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
			final List<String> args= new ArrayList<>();
			args.add("-R"); args.add(support.resource(ref));
			args.add("-n"); args.add("SKIP");
			args.add("-o"); args.add(out.toString());
			if(plain) {
				args.add("--plain");
				}
			args.add(tmp.toString());
			
			
			Assert.assertEquals(
				new VcfSetSequenceDictionary().instanceMain(args)
				,0);
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
