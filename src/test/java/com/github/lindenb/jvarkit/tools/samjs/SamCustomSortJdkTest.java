package com.github.lindenb.jvarkit.tools.samjs;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class SamCustomSortJdkTest extends TestUtils {
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new ParamCombiner().
			initList(collectAllSamOrBam()).
			product(
					"return R1.getMappingQuality() - R2.getMappingQuality();"

				).
			build();
		}
	@Test(dataProvider="src1")
	public void test1(final String inBam,final String expr) throws IOException {
		final File out = createTmpFile(".bam");
		Assert.assertEquals(0,new SamCustomSortJdk().instanceMain(newCmd().add(
        		"-o",out.getPath(),
        		"-e",expr,
        		inBam
        		).
				make()));
		assertIsValidBam(out);
		}
	}
