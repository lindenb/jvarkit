package com.github.lindenb.jvarkit.util.vcf;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.function.Predicate;
import java.util.stream.Stream;

import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class JexlVariantPredicateTest {
	private final TestSupport support = new TestSupport();

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		Stream<String> L1 = support.allVcfOrBcf();
		Stream<String> L2 = Arrays.asList("",
				"QUAL > 0",
				"CHROM == \"chr1\""
				).stream();
		
		return support.combine2(L1, L2);
	}

@Test(dataProvider="src1")
public void test01(final String inputFile,final String expr) 
	throws IOException
	{
	Predicate<VariantContext> pred = JexlVariantPredicate.create(expr);
	final VCFFileReader r =new VCFFileReader(Paths.get(inputFile),false);
	r.iterator().stream().filter(pred);
	r.close();
	}

}
