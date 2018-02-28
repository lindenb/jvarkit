package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.DataProvider;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class Biostar78285Test extends TestUtils{
	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][] {
			{"./src/test/resources/toy.bam"},
			{"--bed /src/test/resources/toy.bed.gz -m 5 -m 10 -m 15 ./src/test/resources/toy.bam"},
			{"--bed /src/test/resources/toy.bed.gz --gcp 10 -R ./src/test/resources/toy.fa -m 5 -m 10 -m 15 ./src/test/resources/toy.bam"}
			};
		}
	
	public void test1(final String cmdsrc) throws IOException {
		final File out = createTmpFile(".vcf");
		Assert.assertEquals(
			new Biostar78285().instanceMain(newCmd().
			add("-o").add(out).
			split(cmdsrc).
			make()
			),0);
		assertIsVcf(out);
		}
}
