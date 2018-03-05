package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.variant.vcf.VCFFileReader;

public class VcfCutSamplesTest extends TestUtils {

@Test(dataProvider="all-vcf-files")
public void test01(final String vcfpath) throws IOException {
	VCFFileReader r= new VCFFileReader(new File(vcfpath),false);
	final List<String> samples = new ArrayList<>(r.getFileHeader().getSampleNamesInOrder());
	r.close();	
	if(samples.isEmpty() ) return;
	Collections.shuffle(samples, this.random);
	final List<String> keep =  samples.subList(0,
			super.random.nextInt(samples.size())
			);
	
	final File vcfOut =  super.createTmpFile(".vcf");
	Assert.assertEquals(new VcfCutSamples().instanceMain(newCmd().
		add("-o",vcfOut.getPath()).
		split(keep.stream().map(S->"-S "+S).collect(Collectors.joining(" "))).
		add(vcfpath).make()
		),0);
	assertIsVcf(vcfOut);
	
	r= new VCFFileReader(vcfOut,false);
	Assert.assertEquals(new TreeSet<>(keep), new TreeSet<>(r.getFileHeader().getSampleNamesInOrder()));
	r.close();
	}
}
