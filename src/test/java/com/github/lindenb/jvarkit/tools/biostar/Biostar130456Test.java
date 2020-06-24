package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.variant.vcf.VCFReader;

public class Biostar130456Test {
	private final TestSupport support = new TestSupport();

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return (Object[][])support.
				allVcfOrBcf().
				map(F->new Object[] {F}).
				toArray()
				;
		}
	
@Test(dataProvider="src1")
public void test01(final String vcfpath) throws IOException {
	try {
	final VCFReader r= VCFReaderFactory.makeDefault().open(new File(vcfpath),false);
	final Set<String> samples = new HashSet<>(r.getHeader().getSampleNamesInOrder());
	r.close();	
	if(samples.isEmpty()) return;
	
	final String vcfOut = IOUtils.getDefaultTmpDir().getPath() + 
			File.separatorChar+
			"tmp.__SAMPLE__.vcf.gz";
	final Path tsvout = support.createTmpPath(".txt");
	Assert.assertEquals(new Biostar130456().instanceMain(new String[] {
		"-o",tsvout.toString(),
		"-p",vcfOut,
		vcfpath
		}),0);
	Assert.assertTrue(Files.exists(tsvout));
	for(final String s:samples) {
		final Path vcfIn2 = Paths.get(vcfOut.replaceAll("__SAMPLE__", s));
		support.assertIsVcf(vcfIn2);
		final VCFReader r2= VCFReaderFactory.makeDefault().open(vcfIn2,false);
		final List<String> samples2 = r2.getHeader().getSampleNamesInOrder();
		r2.close();
		Assert.assertTrue(samples2.size()==1);
		Assert.assertEquals(samples2.get(0), s);
		Assert.assertTrue(Files.deleteIfExists(vcfIn2));
		}
	} finally {
		support.removeTmpFiles();
		}
	}
	
}
