package com.github.lindenb.jvarkit.tools.tview;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;


@AlsoTest(LauncherTest.class)
public class TViewCmdTest {

	private final TestSupport support = new TestSupport();

	
	@DataProvider(name="rf_regions")
public Object[][] getDataRegions() throws IOException  {
	return support.toArrayArray(
			support.randomIntervalsFromDict(Paths.get(support.resource("rotavirus_rf.dict")),10,1000).
			stream().
			map(I->I.getContig()+":"+I.getStart()+"-"+I.getEnd()).
			map(S->new Object[]{S})
			)
			;
	}
	
@Test(dataProvider="rf_regions")
public void test01(final String rgn) throws IOException {
		try {
		final Path imgOut = support.createTmpPath(".txt");
		List<String> args= new ArrayList<>();
		args.add("-V");
		args.add(support.resource("rotavirus_rf.vcf.gz"));
		args.add("-R");
		args.add(support.resource("rotavirus_rf.vcf.fa"));
		args.add("-r");
		args.add(rgn);
		args.add("-o");
		args.add(imgOut.toString());
		for(int i=1;i<=5;i++) {
			args.add(support.resource("S"+i+".bam"));
		}
		
		Assert.assertEquals(new  TViewCmd().instanceMain(args),0);
		support.assertIsNotEmpty(imgOut);
		} finally {
		support.removeTmpFiles();
	}
	}
}
