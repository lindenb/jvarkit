package com.github.lindenb.jvarkit.tools.gnomad;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.tabix.TabixIndex;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

@AlsoTest(LauncherTest.class)
public class VcfGnomadTest {
	private final TestSupport support = new TestSupport();

	
@DataProvider(name="src01")
public Object[][] testData01() {
	return new Object[][] {
		{support.resource("test_vcf01.vcf")}
		
	};
}

private Path createManifest() throws IOException {
	final Path manifestFile = support.createTmpPath(".mft");
	final PrintWriter pw = new PrintWriter(Files.newOutputStream(manifestFile));
	pw.println("exome\t*\t"+support.resource("gnomad.exomes.r2.0.1.sites.vcf.gz"));
	pw.println("genome\t*\t"+support.resource("gnomad.genomes.r2.0.1.sites.1.vcf.gz"));
	pw.flush();
	pw.close();
	support.assertTsvTableIsConsitent(manifestFile, null);
	return manifestFile;
	}

@Test(dataProvider="src01")
public void test01(final String vcfpath) throws IOException {
	try {
		final Path mFile = createManifest();
		final Path vcfOut = support.createTmpPath(".vcf");
		
		Assert.assertEquals(new VcfGnomad().instanceMain(new String[]{
				"-o",vcfOut.toString(),
				"-m",mFile.toString(),
				"--gnomadFilter","MYF111",
				vcfpath
				}),0);
		support.assertIsVcf(vcfOut);
		}
	finally
		{
		support.removeTmpFiles();
		}
	}

@Test(dataProvider="src01")
public void testHeaderV2_1(final String vcfpath) throws IOException {
	try {
		
		final Path manifestFile = support.createTmpPath(".mft");
		final PrintWriter pw = new PrintWriter(Files.newOutputStream(manifestFile));
		final SAMSequenceDictionary dict= SequenceDictionaryUtils.extractRequired(Paths.get(support.resource("human_b37.dict")));

		for(int i=0;i< 2;++i)
			{
			final Path gf = support.createTmpPath(".vcf.gz");
			final PrintWriter pw2 = new PrintWriter(new BlockCompressedOutputStream(Files.newOutputStream(gf),gf));
			pw2.println("##fileformat=VCFv4.2");
			pw2.println("##INFO=<ID=AC_nfe,Number=A,Type=Integer,Description=\"\">");
			pw2.println("##INFO=<ID=AN_nfe,Number=A,Type=Integer,Description=\"\">");
			pw2.println("##INFO=<ID=AF_nfe,Number=A,Type=Float,Description=\"\">");
			pw2.println("##INFO=<ID=non_neuro_AC_nfe_male,Number=A,Type=Integer,Description=\"\">");
			pw2.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
			pw2.flush();
			pw2.close();
			
			support.assertIsVcf(gf);
			
			final TabixIndex idx = IndexFactory.createTabixIndex(gf.toFile(), new VCFCodec(),dict);
			Path gfidx= Paths.get(gf.toString()+".tbi");
			
			support.deleteOnExit(gfidx);
			idx.write(gfidx);
			
			pw.println(""+ (i==0?"genome":"exome")+"\t*\t"+ gf.toString());
			}
		pw.flush();
		pw.close();
		
		final Path vcfOut = support.createTmpPath(".vcf");
		
		Assert.assertEquals(new VcfGnomad().instanceMain(new String[]{
				"-o",vcfOut.toString(),
				"--ani",
				"-m",manifestFile.toString(),
				vcfpath
				}),0);
		support.assertIsVcf(vcfOut);
		try (VCFFileReader vcg=new VCFFileReader(vcfOut, false)) {
			VCFHeader h = vcg.getFileHeader();
			Assert.assertNotNull(h);
			for(int i=0;i< 2;++i) {
				String prefix="gnomad_"+(i==0?"exome":"genome")+"_";
				VCFInfoHeaderLine ih = h.getInfoHeaderLine(prefix+"AF_NFE");
				Assert.assertNotNull(ih);
				Assert.assertEquals(ih.getCountType(),VCFHeaderLineCount.A);
				Assert.assertEquals(ih.getType(),VCFHeaderLineType.Float);

	
				ih = h.getInfoHeaderLine(prefix+"AC_NFE");
				Assert.assertNotNull(ih);
				Assert.assertEquals(ih.getCountType(),VCFHeaderLineCount.A,"Type is "+ih.getCountType());
				Assert.assertEquals(ih.getType(),VCFHeaderLineType.Integer);

				
				ih = h.getInfoHeaderLine(prefix+"AN_NFE");
				Assert.assertNotNull(ih);
				Assert.assertTrue(ih.isFixedCount(),ih.toString());
				Assert.assertEquals(ih.getType(),VCFHeaderLineType.Integer,"Type is "+ih.getCountType());
				Assert.assertEquals(ih.getCount(),1);
				}
			}
		}
	finally
		{
		support.removeTmpFiles();
		}
	}


}
