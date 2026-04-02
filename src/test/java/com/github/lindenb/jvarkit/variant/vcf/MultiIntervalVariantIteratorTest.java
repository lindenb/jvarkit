package com.github.lindenb.jvarkit.variant.vcf;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class MultiIntervalVariantIteratorTest {
	@Test
	public void testMultiIterator() throws Exception {
		final TestSupport support=new TestSupport();
		try {
			Path vcf = Paths.get(support.resource("rotavirus_rf.vcf.gz"));
			Assert.assertTrue(Files.exists(vcf));
			try(VCFFileReader vcfFileReader = new VCFFileReader(vcf, true)) {
				try(CloseableIterator<VariantContext> iter= MultiIntervalVariantIterator.query(vcfFileReader, Arrays.asList(
						new SimpleInterval("RF01", 1, 10000),
						new SimpleInterval("_no_in_dict", 1, 10000),
						new SimpleInterval("RF01", 1, 10000),
						new SimpleInterval("RF02", 1, 10000),
						new SimpleInterval("RF03", 1, 10000),
						new SimpleInterval("RF06", 1129, 1129),
						new SimpleInterval("RF06", 1131, 1131)
						))) {
					Assert.assertTrue(iter.hasNext());
					final List<VariantContext> L= iter.toList();
					Assert.assertTrue(L.stream().anyMatch(R->R.overlaps(new SimpleInterval("RF06", 1129, 1129))));
					Assert.assertFalse(L.isEmpty());
					for(int i=0;i+1< L.size();i++) {
						final VariantContext vc1 = L.get(i);
						for(int j=i+1; j< L.size();j++) {
							final VariantContext vc2 = L.get(j);
							Assert.assertFalse(
								vc1.contigsMatch(vc2) &&
								vc1.getStart()==vc2.getStart() &&
								vc1.getEnd()==vc2.getEnd() &&
								vc1.getAlleles().equals(vc2.getAlleles())
								);
							}
						}
					}
				}
			}
		catch(final Throwable err) {
			Assert.fail("boum",err);
			}
		}
}