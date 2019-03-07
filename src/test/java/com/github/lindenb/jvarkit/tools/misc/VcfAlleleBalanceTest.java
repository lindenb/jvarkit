package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.vcf.VCFUtilsTest;

@AlsoTest(VCFUtilsTest.class)
public class VcfAlleleBalanceTest {
	private final TestSupport support = new TestSupport();

	@DataProvider(name="src01")
	public Object[][] testData01() {
			return support.toArrayArray(
					support.allVcfOrBcf().
					map(S->new Object[] {S})
					);

			}
	
	@Test(dataProvider="src01")
	public void test01(final String inputFile) 
		throws IOException
		{
		try {
		final Path out = support.createTmpPath(".vcf");
		Assert.assertEquals(0,new VcfAlleleBalance().instanceMain(new String[] {
			"-o",out.toString(),
			inputFile
			}));
		support.assertIsVcf(out);
		} finally {
			support.removeTmpFiles();
		}
		}
	@Test(dataProvider="src01")
	public void testWithPed(final String inputFile) 
		throws IOException
		{
		try {
		final Path out = support.createTmpPath(".vcf");
		final Path ped = support.createRandomPedigreeFromFile(inputFile);
		if(ped==null) return;
		Assert.assertEquals(0,new VcfAlleleBalance().instanceMain(new String[] {
			"-o",out.toString(),
			"-p",ped.toString(),
			inputFile
			}));
		support.assertIsVcf(out);
		} finally {
			support.removeTmpFiles();
		}
		}
	}
