package com.github.lindenb.jvarkit.tools.gtf;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReaderTest;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

import htsjdk.samtools.util.IOUtil;

@AlsoTest({LauncherTest.class,GtfReaderTest.class})
public class GtfFileSplitterTest {
	private final TestSupport support = new TestSupport();
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		final String gtf1= support.resource("gencode.v19.annotation.gtf");
		final String gtf2= support.resource("Homo_sapiens.GRCh37.87.gtf.gz");
		return  new Object[][] {
				{"gene",gtf1},
				{"transcript",gtf1},
				{"contig",gtf1},
				{"group",gtf1},
				{"stack",gtf1},
				
				{"gene",gtf2},
				{"transcript",gtf2},
				{"contig",gtf2},
				{"group",gtf2},
				{"stack",gtf2}
				};
		}
	
	@Test(dataProvider="src1")
	public void runMethod(final String method,String gtf) throws IOException {
		File tmp=null;
		try {
			tmp = IOUtil.createTempDir("tmp.", ".dir");
			
			Assert.assertEquals(new GtfFileSplitter().instanceMain(new String[] {
					"-o",tmp.toString(),
					"-m",method,
					gtf
					}),0);
			
			}
		finally {
			if(tmp!=null) IOUtil.deleteDirectoryTree(tmp);
			support.removeTmpFiles();
			}
	}
	
}
