package com.github.lindenb.jvarkit.tools.biostar;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class Biostar139647Test {
	private final TestSupport support = new TestSupport();
	@Test
	public void test01() throws IOException {
		try {
		final Path datain = support.createTmpPath(".fa");
		final PrintWriter pw =new PrintWriter(Files.newOutputStream(datain));
		pw.println(">A\n--ATGACT-ATACTAC-");
		pw.println(">B\nAAACTGACAATA-TACA");
		pw.println(">C\n--ATGACT-ATACTAC-");
		pw.flush();
		pw.close();
		
		final Path out = support.createTmpPath(".bam");
		Assert.assertEquals(new Biostar139647().instanceMain(new String[] {
			"-o",out.toString(),
			datain.toString()
			}),0);
		support.assertIsValidBam(out);
		}
	finally
		{
		this.support.removeTmpFiles();
		}
	}
}
