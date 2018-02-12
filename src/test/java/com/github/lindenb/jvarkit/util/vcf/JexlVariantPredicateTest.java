package com.github.lindenb.jvarkit.util.vcf;

import java.io.File;
import java.io.IOException;
import java.util.function.Predicate;

import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class JexlVariantPredicateTest extends TestUtils{

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new ParamCombiner().
			initList(collectAllVcfs()).
			product(
					"",
					"QUAL > 0",
					"CHROM == \"chr1\""
					).
			build();
	}

@Test(dataProvider="src1")
public void test01(final String inputFile,final String expr) 
	throws IOException
	{
	Predicate<VariantContext> pred = JexlVariantPredicate.create(expr);
	final VCFFileReader r =new VCFFileReader(new File(inputFile),false);
	r.iterator().stream().filter(pred);
	r.close();
	}

}
