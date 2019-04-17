package com.github.lindenb.jvarkit.tools.vcf2sql;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class VcfToSqlTest {
	
	private final TestSupport support =new TestSupport();

	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return support.toArrayArray(support.
				allVcfOrBcf().
				map(F->new Object[] {F})
				)
				;
		}
	
	@Test(dataProvider = "src1")
	public void test01(final String vcf) throws IOException {
		try {
			final Path sqlout = support.createTmpPath(".sql");
			Assert.assertEquals(new VcfToSql().instanceMain(new String[] {
				"-o",sqlout.toString(),
				vcf
				}),0);
			 support.assertIsNotEmpty(sqlout);
		} finally {
			support.removeTmpFiles();
		}
	}
}
