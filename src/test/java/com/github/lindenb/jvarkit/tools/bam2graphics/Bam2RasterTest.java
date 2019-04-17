package com.github.lindenb.jvarkit.tools.bam2graphics;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class Bam2RasterTest {
	private final TestSupport support = new TestSupport();

	@DataProvider(name="rf_regions")
	public Object[][] getDataRegions() throws IOException  {
		return new Object[][] {
			{"RF01:1-100"},
			{"RF02:1-100"},
			{"RF03:1-100"}
			}
			;
		}

private void run(final String rgn,Path out) throws IOException{
	final List<String> args= new ArrayList<>();
	args.add("-R");
	args.add(support.resource("rotavirus_rf.fa"));
	args.add("-r");
	args.add(rgn);
	args.add("-o");
	args.add(out.toString());
	for(int i=1;i<=5;i++) args.add(support.resource("S"+i+".bam"));
	Assert.assertEquals(new Bam2Raster().instanceMain(args),0);
	}
	
@Test(dataProvider="rf_regions")
public void test01(final String rgn) throws IOException {
	try {
		final Path imgOut = support.createTmpPath(".png");
		run(rgn,imgOut);
		support.assertIsNotEmpty(imgOut);
		}
	finally
		{
		support.removeTmpFiles();
		}
	}

@Test(dataProvider="rf_regions")
public void testZip(final String rgn) throws IOException {
	try {
		final Path imgOut = support.createTmpPath(".zip");
		run(rgn,imgOut);
		support.assertZip(imgOut);
		}
	finally
		{
		support.removeTmpFiles();
		}
	}

}
