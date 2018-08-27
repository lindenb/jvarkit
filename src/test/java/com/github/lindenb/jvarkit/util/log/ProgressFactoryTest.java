package com.github.lindenb.jvarkit.util.log;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
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
@Test(dataProvider="all-sam-or-bam-files")
void testBam01(final String bam) throws IOException
	{
	final SamReader sr =  SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT).open(SamInputResource.of(bam));
	Assert.assertTrue(ProgressFactory.newInstance().stream(sr.iterator()).count()>0);
	sr.close();
	}

}
