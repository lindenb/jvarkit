package com.github.lindenb.jvarkit.variant.vcf;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.variant.variantcontext.VariantContext;

public class VCFByIndexTest {
	
	private void testByIndex(String fname,String indexes) throws Exception {
		final TestSupport support=new TestSupport();
		try {
			Path vcf = Paths.get(support.resource(fname));
			Assert.assertTrue(Files.exists(vcf));
			Path ith = support.createTmpPath(VCFByIndex.INDEX_SUFFIX);
			
			VCFByIndex.buildIndex(vcf,ith);
			Assert.assertTrue(VCFByIndex.indexExists(ith.toString()));
			Assert.assertTrue(Files.exists(ith));
			Assert.assertTrue(Files.size(ith)>0L);
			Assert.assertEquals(Files.size(ith)%Long.BYTES,0);
			
			long n_variants = support.wc(vcf);
			Assert.assertTrue(n_variants>0L);
			
			try(VCFByIndex byidx = new VCFByIndex(vcf, ith)) {
				Assert.assertNotNull(byidx.getHeader());
				Assert.assertEquals(byidx.size(), n_variants);
				for(String s: CharSplitter.COMMA.split(indexes)) {
					final long offset = Long.parseLong(s);
					final VariantContext ctx = byidx.get(offset);
					Assert.assertNotNull(ctx);
					}
				}
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
		testByIndex("test_vcf01.vcf","1,10,50,70,90");
		}
	@Test
	public void testVcfGz() throws Exception {
		testByIndex("rotavirus_rf.vcf.gz","10,20,30,40");
		}
	

}
