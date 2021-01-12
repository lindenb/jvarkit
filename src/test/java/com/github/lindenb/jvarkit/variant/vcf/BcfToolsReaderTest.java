package com.github.lindenb.jvarkit.variant.vcf;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;

@AlsoTest(BcfToolsUtilsTest.class)
public class BcfToolsReaderTest {
private final TestSupport support  = new TestSupport();
@Test
public void test1() throws IOException {
	try {
		final Path p = Paths.get(support.resource("toy.bcf"));
		Assert.assertTrue(BcfToolsUtils.isBcfToolsRequired(p));
		int n=0;
		if(BcfToolsUtils.isBcfToolsRequired(p)) {
			try(BcfToolsReader br = new BcfToolsReader(p.toString())) {
				br.getHeader();
				try(CloseableIterator<VariantContext> iter=br.iterator()) {
					while(iter.hasNext()) {
						iter.next();
						n++;
						}
					}
				}
			Assert.assertEquals(n, 5);
			}
		}
	finally
		{
		support.removeTmpFiles();
		}
	}
}
