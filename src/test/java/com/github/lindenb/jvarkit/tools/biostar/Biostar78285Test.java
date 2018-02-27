package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class Biostar78285Test extends TestUtils{
	
	public void test1() throws IOException {
		final File out = createTmpFile(".vcf");
		Assert.assertEquals(
			new Biostar78285().instanceMain(newCmd().
			add("-o").add(out).
			add("--bed","/src/test/resources/toy.bed.gz").
			add("-m").add(5).
			add("-m").add(10).
			add("-m").add(15).
			add("./src/test/resources/toy.bam").
			make()
			),0);
		assertIsVcf(out);
		}
}
