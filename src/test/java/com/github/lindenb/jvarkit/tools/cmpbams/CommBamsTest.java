package com.github.lindenb.jvarkit.tools.cmpbams;

import java.io.File;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;


public class CommBamsTest extends TestUtils {

@Test(dataProvider="all-sam-or-bam-files")
public void test1(final String bampath) throws Exception
	{
	final File out = super.createTmpFile(".txt");
	final File sortedBam1= sortBamOnQueryName(Paths.get(bampath),R->random.nextDouble()<0.5);
	final File sortedBam2= sortBamOnQueryName(Paths.get(bampath),R->random.nextDouble()<0.5);
	Assert.assertEquals(new CommBams().instanceMain(new String[] {
		"-o",out.getPath(),
		sortedBam1.getPath(),
		sortedBam2.getPath(),
		}),0);	
	}
}
