package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.variant.vcf.VCFFileReader;

public class Biostar130456Test extends TestUtils{

@Test(dataProvider="all-vcf-files")
public void test01(final String vcfpath) throws IOException {
	final VCFFileReader r= new VCFFileReader(new File(vcfpath),false);
	final Set<String> samples = new HashSet<>(r.getFileHeader().getSampleNamesInOrder());
	r.close();	
	if(samples.isEmpty()) return;
	
	final String vcfOut = IOUtils.getDefaultTmpDir().getPath() + 
			File.separatorChar+
			"tmp.__SAMPLE__.vcf.gz";
	final File tsvout = createTmpFile(".txt");
	Assert.assertEquals(new Biostar130456().instanceMain(new String[] {
		"-o",tsvout.getPath(),
		"-p",vcfOut,
		vcfpath
		}),0);
	Assert.assertTrue(tsvout.exists());
	for(final String s:samples) {
		final File vcfIn2 = new File(vcfOut.replaceAll("__SAMPLE__", s));
		assertIsVcf(vcfIn2);
		final VCFFileReader r2= new VCFFileReader(vcfIn2,false);
		final List<String> samples2 = r2.getFileHeader().getSampleNamesInOrder();
		r2.close();
		Assert.assertTrue(samples2.size()==1);
		Assert.assertEquals(samples2.get(0), s);
		Assert.assertTrue(vcfIn2.delete());
		}
	}
	
}
