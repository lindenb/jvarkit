package com.github.lindenb.jvarkit.tools.vcf2sql;

import java.io.File;
import java.io.IOException;
import org.testng.Assert;
import org.testng.annotations.Test;
import com.github.lindenb.jvarkit.tools.tests.TestUtils;


public class VcfToSqlTest extends TestUtils {
	@Test(dataProvider = "all-vcf-files")
	public void test01(final String vcf) throws IOException {
		final File sqlout = createTmpFile(".sql");
		Assert.assertEquals(new VcfToSql().instanceMain(new String[] {
			"-o",sqlout.getPath(),
			vcf
			}),0);
		assertIsNotEmpty(sqlout);
	}
}
