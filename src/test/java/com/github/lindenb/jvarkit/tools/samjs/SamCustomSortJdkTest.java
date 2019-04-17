package com.github.lindenb.jvarkit.tools.samjs;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class SamCustomSortJdkTest  {
	
	private final TestSupport support = new TestSupport();

	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		List<String> L1= support.allSamOrBams().collect(Collectors.toList());
		List<String> L2= Collections.singletonList("return R1.getMappingQuality() - R2.getMappingQuality();");
		return support.combine2(L1.stream(), L2.stream());
		}
	@Test(dataProvider="src1")
	public void test1(final String inBam,final String expr) throws IOException {
		try {
			final Path out = support.createTmpPath(".bam");
			Assert.assertEquals(0,new SamCustomSortJdk().instanceMain(new String[] {
        		"-o",out.toString(),
        		"-e",expr,
        		inBam
				}));
			support.assertIsValidBam(out);
			}
		finally
			{
			support.removeTmpFiles();
			}
		}
	}
