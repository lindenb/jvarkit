package com.github.lindenb.jvarkit.tools.biostar;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class Biostar95652Test{
	private final TestSupport support = new TestSupport();
	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][]{
			{"NP_077719.2  XP_513697.3 XP_001114248.1 XP_540266.3 XP_002686160.2 NP_035058.2 NP_077334.1  NP_001238962.1 NP_001108566.1"}
			};
		}
	
	@Test(dataProvider="src1")
	public void test1(final String acns) throws IOException {
		try {
			final Path out = support.createTmpPath(".svg");
			final List<String> args = new ArrayList<>();
			args.add("-o");
			args.add(out.toString());
			args.addAll(Arrays.asList(acns.split("[ ]+")));
			
			Assert.assertEquals( new Biostar95652().instanceMain(args),0);
			support.assertIsXml(out);
			}
		finally {
			support.removeTmpFiles();
			}
		}

}
