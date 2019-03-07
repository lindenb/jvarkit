package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;

public class SamReadLengthDistributionTest  {
@DataProvider(name="src01")
public Object[][] src01(){
	final SAMRecordPartition v[]=SAMRecordPartition.values();
	final Object[][] ret=new Object[v.length][];
	for(int i=0;i< v.length;i++) ret[i]=new Object[] {v[i].name()};
	return ret;
}
	
@Test(dataProvider="src01")
public void test01(final String partition) throws IOException{
	final TestSupport support = new TestSupport();
	try {
		final Path out = support.createTmpPath(".txt");
		
		Assert.assertEquals(
			new SamReadLengthDistribution().instanceMain(new String[] {
			"-o",out.toString(),
			"--groupby",partition,
			support.resource("S1.bam"),
			support.resource("S2.bam"),
			support.resource("S3.bam"),
			support.resource("S4.bam"),
			support.resource("S5.bam")
			}),0);
		support.assertTsvTableIsConsitent(out, null);
		}
	finally {
		support.removeTmpFiles();
		}
	}
}
