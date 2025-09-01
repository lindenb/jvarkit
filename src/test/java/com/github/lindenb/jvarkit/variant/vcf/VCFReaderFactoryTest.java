package com.github.lindenb.jvarkit.variant.vcf;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Iterator;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFReader;

public class VCFReaderFactoryTest {
@DataProvider(name = "src1")
public Iterator<Object[]> listVCFsFortest() {
	return Arrays.asList(
			"toy.bcf",
			"toy.vcf.gz"
			).stream()
		.map(it->new Object[] {it})
		.iterator();
	}

@Test(dataProvider = "src1")
public void scanFile(String fname) throws IOException {
	final TestSupport support  = new TestSupport();
	long n=0;
	final Path p = Paths.get(support.resource(fname));
	try(VCFReader r =VCFReaderFactory.makeDefault().open(p,!fname.endsWith(FileExtensions.VCF))) {
		r.getHeader();
		try(CloseableIterator<VariantContext> iter=r.iterator()) {
			while(iter.hasNext()) {
				iter.next();
				n++;
				}
			}
		}
	Assert.assertEquals(n, 5);
	}
		

}
