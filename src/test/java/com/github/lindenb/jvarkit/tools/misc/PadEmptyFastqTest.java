package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class PadEmptyFastqTest {
	private final TestSupport support = new TestSupport();

	
	@DataProvider(name="src01")
	public Object[][] getData() {
		List<Object[]> L=new ArrayList<>();
		for(int s=1;s<=5;s++) {
			for(int t=1;t<=2;t++) {
				L.add(new Object[] {support.resource("S"+s+".R"+(t)+".fq.gz")});
			}
		}
		
		return support.toArrayArray(L.stream());
	}
	
	@Test(dataProvider="src01")
	public void test01(final String fq) 
			throws IOException
			{
			try {
				final Path out = support.createTmpPath(".fq");
			
				Assert.assertEquals(0,new PadEmptyFastq().instanceMain(new String[] {
					"-o",out.toString(),
					fq
					}));
				support.assertIsFastq(out);
				}
			finally {
				support.removeTmpFiles();
			}
			}
}
