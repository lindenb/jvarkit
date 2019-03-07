package com.github.lindenb.jvarkit.tools.vcfsparql;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.vcf.VCFUtilsTest;


@AlsoTest(VCFUtilsTest.class)
public class VcfSparqlTest {
	private final TestSupport support = new TestSupport();

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		List<String> L1= Arrays.asList(support.resource("S1.vcf.gz"),support.resource("rotavirus_rf.vcf.gz"));
		List<String> L2= Arrays.asList(
				"SELECT distinct *\n WHERE {  ?s ?p ?o .\n  }",
				"select ?a ?b ?c  {?a a <"+VariantGraph.NS+"Genotype> . ?a ?b ?c .}"
				);
		List<Object[]> L3 = new ArrayList<>();
		for(String s1:L1)
			for(String s2:L2)
				L3.add(new Object[] {s1,s2});
		
		return support.toArrayArray(L3.stream());
		}
	
	@Test(dataProvider="src1")
	public void test1(final String vcf,final String query) throws IOException {
			try {
			final Path out = support.createTmpPath(".txt");
			Assert.assertEquals(
				new VcfSparql().instanceMain(new String[] {
				"-a",
				"-o",out.toString(),
				"-r","RF01:1-1000",
				"-e",query,
				vcf}
				),0);
			support.assertIsNotEmpty(out);
			}
		finally {
			support.removeTmpFiles();
			}
		}

}
