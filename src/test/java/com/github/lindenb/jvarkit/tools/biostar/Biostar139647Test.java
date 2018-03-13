package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;


public class Biostar139647Test extends TestUtils {
	@Test
	public void test01() throws IOException {
		
		final File datain = createTmpFile(".fa");
		final PrintWriter pw =new PrintWriter(datain);
		pw.println(">A\n--ATGACT-ATACTAC-");
		pw.println(">B\nAAACTGACAATA-TACA");
		pw.println(">C\n--ATGACT-ATACTAC-");
		pw.flush();
		pw.close();
		
		final File out = createTmpFile(".bam");
		Assert.assertEquals(new Biostar139647().instanceMain(new String[] {
			"-o",out.getPath(),
			datain.getPath()
			}),0);
		assertIsValidBam(out);
	}
}
