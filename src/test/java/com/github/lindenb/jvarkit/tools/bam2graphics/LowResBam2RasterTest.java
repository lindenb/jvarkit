package com.github.lindenb.jvarkit.tools.bam2graphics;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class LowResBam2RasterTest {
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
	
@Test(dataProvider="rf_regions")
public void test01(final String rgn) throws IOException {
	try {
		Path imgOut = support.createTmpPath(".png");
		
		final List<String> args= new ArrayList<>();
		args.add("-kg");
		args.add(support.resource("rotavirus_rf.knowngenes.tsv.gz"));
		args.add("-R");
		args.add(support.resource("rotavirus_rf.fa"));
		args.add("-r");
		args.add(rgn);
		args.add("-o");
		args.add(imgOut.toString());
		Arrays.asList("1","2","3","4","5").stream().
		map(S->support.resource("/S"+S+".bam")).
		forEach(S->args.add(S));
		
		Assert.assertEquals(new LowResBam2Raster().instanceMain(
				args.toArray(new String[args.size()])),
				0);
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
		final List<String> args= new ArrayList<>();
		args.add("-kg");
		args.add(support.resource("rotavirus_rf.knowngenes.tsv.gz"));
		args.add("-R");
		args.add(support.resource("rotavirus_rf.fa"));
		args.add("-o");
		args.add(imgOut.toString());
		Arrays.asList("1","2","3","4","5").stream().
			map(S->support.resource("/S"+S+".bam")).
			forEach(S->args.add(S));
		

		Assert.assertEquals(new LowResBam2Raster().instanceMain(
			args.toArray(new String[args.size()])),
			0);
		support.assertZip(imgOut);
		}
	finally
		{
		support.removeTmpFiles();
		}
	}

}
