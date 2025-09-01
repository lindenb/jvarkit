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

import htsjdk.variant.vcf.VCFIterator;

public class BcfIteratorBuilderTest {
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
public void scanIterator(String fname) throws IOException {
	final TestSupport support  = new TestSupport();
	long n=0;
	final Path p = Paths.get(support.resource(fname));
	try(VCFIterator iter = new BcfIteratorBuilder().open(p.toString())) {
		iter.getHeader();
		while(iter.hasNext()) {
			iter.next();
			n++;
			}
		}
	Assert.assertEquals(n, 5);
	}
		
	
	
}
