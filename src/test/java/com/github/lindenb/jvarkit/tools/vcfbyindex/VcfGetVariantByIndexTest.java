package com.github.lindenb.jvarkit.tools.vcfbyindex;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class VcfGetVariantByIndexTest {
	
	private void testByIndex(String fname,boolean compressed,String indexes) throws Exception {
		final TestSupport support=new TestSupport();
		try {
			Path vcf1= support.createTmpPath(compressed?".vcf.gz":".vcf");
			Path outvcf= support.createTmpPath(".vcf");
			Path invcf = Paths.get(support.resource(fname));
			Assert.assertTrue(Files.exists(invcf));
			IOUtils.copyTo(invcf, vcf1);
			Assert.assertEquals(new VcfGetVariantByIndex().instanceMain(new String[] {
					"--force",
					"--index",indexes,
					"-o",outvcf.toString(),
					vcf1.toString()}),0);
			final Path ith = vcf1.getParent().resolve(vcf1.getFileName().toString()+".ith");
			Assert.assertTrue(Files.exists(ith),"index of "+vcf1+" is "+ith+" (should exist)");
			support.deleteOnExit(ith);
			support.assertIsVcf(outvcf);
			Assert.assertEquals(support.wc(outvcf),CharSplitter.COMMA.split(indexes).length,"wc("+outvcf+") vs "+indexes);
			}
		catch(Throwable err) {
			Assert.fail(fname, err);
			}
		finally {
			support.removeTmpFiles();
			}
		}
	@Test
	public void testVcfPlain() throws Exception {
		testByIndex("test_vcf01.vcf",false,"1,10,50,70,90");
		}
	@Test
	public void testVcfGz() throws Exception {
		testByIndex("rotavirus_rf.vcf.gz",true,"10,20,30,40");
		}
	}
