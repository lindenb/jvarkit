package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class Biostar173114Test extends TestUtils{
private  File basetest(final String bam,String param) throws IOException {
	final File out = createTmpFile(".bam");
	Assert.assertEquals(
		new Biostar173114().instanceMain(newCmd().
		add("-o").add(out).
		split(param).
		add(bam).
		make()
		),0);
	assertIsValidBam(out);
	return out;
	}

@Test(dataProvider="all-sam-or-bam-files")
public void test01(final String bam)throws IOException {
	basetest(bam,"");
	}
@Test(dataProvider="all-sam-or-bam-files")
public void keepSequence(final String bam)throws IOException {
	basetest(bam,"--keepSequence");
	}
@Test(dataProvider="all-sam-or-bam-files")
public void keepName(final String bam)throws IOException {
	basetest(bam,"--keepName");
	}
@Test(dataProvider="all-sam-or-bam-files")
public void keepAttributes(final String bam)throws IOException {
	basetest(bam,"--keepAttributes");
	}
@Test(dataProvider="all-sam-or-bam-files")
public void keepReadGroup(final String bam)throws IOException {
	basetest(bam,"--keepReadGroup");
	}
@Test(dataProvider="all-sam-or-bam-files")
public void mate(final String bam)throws IOException {
	basetest(bam,"--mate");
	}
@Test(dataProvider="all-sam-or-bam-files")
public void keepQualities(final String bam)throws IOException {
	basetest(bam,"--keepQualities");
	}
@Test(dataProvider="all-sam-or-bam-files")
public void keepCigar(final String bam)throws IOException {
	basetest(bam,"--keepCigar");
	}
}
