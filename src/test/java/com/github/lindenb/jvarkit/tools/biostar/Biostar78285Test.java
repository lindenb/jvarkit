package com.github.lindenb.jvarkit.tools.biostar;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.DataProvider;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class Biostar78285Test {
	private final TestSupport support = new TestSupport();
	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][] {
			{support.resource("toy.bam")},
			{"--bed "+support.resource("toy.bed.gz")+" -m 5 -m 10 -m 15 "+support.resource("toy.bam")},
			{"--bed "+support.resource("toy.bed.gz")+" --gcp 10 -R "+ support.resource("toy.bam") +" -m 5 -m 10 -m 15 "+support.resource("toy.bam")}
			};
		}
	
	public void test1(final String cmdsrc) throws IOException {
		try {
			final Path out = support.createTmpPath(".vcf");
			
			List<String> args=new ArrayList<>();
			args.add("-o");
			args.add(out.toString());
			args.addAll(Arrays.asList(cmdsrc.split("[ \t]+")));
			
			Assert.assertEquals(
				new Biostar78285().instanceMain(args),0);
			support.assertIsVcf(out);
			}
		finally
			{
			support.removeTmpFiles();
			}
		}
}
