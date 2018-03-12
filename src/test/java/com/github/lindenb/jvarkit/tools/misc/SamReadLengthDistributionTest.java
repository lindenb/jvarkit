package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;

public class SamReadLengthDistributionTest extends TestUtils {
@DataProvider(name="src01")
public Object[][] src01(){
	final SAMRecordPartition v[]=SAMRecordPartition.values();
	final Object[][] ret=new Object[v.length][];
	for(int i=0;i< v.length;i++) ret[i]=new Object[] {v[i].name()};
	return ret;
}
	
@Test(dataProvider="src01")
public void test01(final String partition) throws IOException{
	final File out = createTmpFile(".txt");
	Assert.assertEquals(
		new SamReadLengthDistribution().instanceMain(newCmd().
		add("-o").add(out).
		add("--groupby").add(partition).
		add(Arrays.asList("1","2","3","4","5").stream().
				map(S->SRC_TEST_RESOURCE+"/S"+S+".bam").
				toArray()).
		make()
		),0);
	assertTsvTableIsConsitent(out, null);
	}
}
