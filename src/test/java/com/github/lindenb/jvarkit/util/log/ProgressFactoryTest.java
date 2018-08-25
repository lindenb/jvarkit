package com.github.lindenb.jvarkit.util.log;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.variant.vcf.VCFFileReader;

public class ProgressFactoryTest extends TestUtils{

@Test
void test01() throws IOException {
	final VCFFileReader r=new VCFFileReader(new File(SRC_TEST_RESOURCE+"/rotavirus_rf.vcf.gz"),false);
	Assert.assertTrue(ProgressFactory.newInstance().dictionary(r).stream(r.iterator()).count()>0);
	r.close();
	}
@Test
void test02() throws IOException {
	final VCFFileReader r=new VCFFileReader(new File(SRC_TEST_RESOURCE+"/rotavirus_rf.vcf.gz"),false);
	Assert.assertTrue(ProgressFactory.newInstance().stream(r.iterator()).count()>0);
	r.close();
	}

}
