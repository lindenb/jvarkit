package com.github.lindenb.jvarkit.tools.biostar;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class Biostar173114Test{
	private final TestSupport support = new TestSupport();

@DataProvider(name = "src1")
public Object[][] createData1() {
	return (Object[][])support.
			allSamOrBams().
			map(F->new Object[] {F}).
			toArray()
			;
	}

	
private  void basetest(final String bam,String param) throws IOException {
	try {
		final Path out = support.createTmpPath(".bam");
		final List<String> args = new ArrayList<>();
		args.add("-o");
		args.add(out.toString());
		if(!param.isEmpty()) args.add(param);
		args.add(bam);
		
		Assert.assertEquals( new Biostar173114().instanceMain(args),0);
		support.assertIsValidBam(out);
		Assert.assertEquals(support.wc(Paths.get(bam)),support.wc(out));
		}
	finally {
		support.removeTmpFiles();
		}
	}

@Test(dataProvider="src1")
public void test01(final String bam)throws IOException {
	basetest(bam,"");
	}
@Test(dataProvider="src1")
public void keepSequence(final String bam)throws IOException {
	basetest(bam,"--keepSequence");
	}
@Test(dataProvider="src1")
public void keepName(final String bam)throws IOException {
	basetest(bam,"--keepName");
	}
@Test(dataProvider="src1")
public void keepAttributes(final String bam)throws IOException {
	basetest(bam,"--keepAttributes");
	}
@Test(dataProvider="src1")
public void keepReadGroup(final String bam)throws IOException {
	basetest(bam,"--keepReadGroup");
	}
@Test(dataProvider="src1")
public void mate(final String bam)throws IOException {
	basetest(bam,"--mate");
	}
@Test(dataProvider="src1")
public void keepQualities(final String bam)throws IOException {
	basetest(bam,"--keepQualities");
	}
@Test(dataProvider="src1")
public void keepCigar(final String bam)throws IOException {
	basetest(bam,"--keepCigar");
	}
}
