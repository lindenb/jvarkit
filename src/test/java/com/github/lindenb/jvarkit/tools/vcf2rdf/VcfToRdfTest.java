package com.github.lindenb.jvarkit.tools.vcf2rdf;

import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class VcfToRdfTest {
	@Test
	public void testVcf2RDF() throws Exception {
		final TestSupport support=new TestSupport();
		try {
			Path ttl= support.createTmpPath(".ttl");
			Assert.assertEquals(new VcfToRdf().instanceMain(new String[] {
					"--hide","",
					"-o",ttl.toString(),
					support.resource("ExAC.r1.sites.vep.vcf.gz")})
					,0);
			
			support.assertIsNotEmpty(ttl);
			}
		catch(Throwable err) {
			Assert.fail("testVcf2RDF", err);
			}
		finally {
			support.removeTmpFiles();
			}
		}
	
}
