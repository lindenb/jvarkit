package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;

import htsjdk.samtools.util.IOUtil;

public class KnownGenesToBedTest extends TestUtils{
	
@DataProvider(name = "src1")
public Object[][] createData1() {
	return new Object[][]{
		{SRC_TEST_RESOURCE+"/rotavirus_rf.knowngenes.tsv.gz"}
		};
	}
@Test(dataProvider="src1")
public void test(final String kgfile) throws IOException {
	File out =super.createTmpFile(".bed");
	Assert.assertEquals(new KnownGenesToBed().instanceMain(new String[] {
			"-o",out.getPath(),
			kgfile
			}),0);
	Assert.assertTrue(out.exists());
	final BedLineCodec codec= new BedLineCodec();
	Assert.assertTrue(IOUtil.slurpLines(out).stream().
			map(L->codec.decode(L)).
			filter(L->L!=null).count()>0L);
	}
}
