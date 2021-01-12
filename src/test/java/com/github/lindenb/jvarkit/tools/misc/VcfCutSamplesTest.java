package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.TreeSet;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;
import com.github.lindenb.jvarkit.util.vcf.VCFUtilsTest;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.variant.vcf.VCFReader;


@AlsoTest({VCFUtilsTest.class,LauncherTest.class})
public class VcfCutSamplesTest  {

	private final TestSupport support = new TestSupport();

	@DataProvider(name="src01")
	public Object[][] testData01() {
			return support.toArrayArray(
					support.allVcfOrBcf().
					map(S->new Object[] {S})
					);
			}
	
	@Test(dataProvider="src01")
	public void test01(final String vcfpath) throws IOException {
		try {
			VCFReader r= VCFReaderFactory.makeDefault().open(Paths.get(vcfpath),false);
			final List<String> samples = new ArrayList<>(r.getHeader().getSampleNamesInOrder());
			r.close();	
			if(samples.isEmpty() ) return;
			Collections.shuffle(samples, support.random);
			final List<String> keep =  samples.subList(0,
					support.random.nextInt(samples.size())
					);
			
			final Path vcfOut =  support.createTmpPath(".vcf");
			List<String> args=new ArrayList<>();
			args.add("-o");
			args.add(vcfOut.toString());
			for(String s : keep) { args.add("-S");args.add(s);}
			args.add(vcfpath);
			
			
			Assert.assertEquals(new VcfCutSamples().instanceMain(args),0);
			support.assertIsVcf(vcfOut);
			
			r= VCFReaderFactory.makeDefault().open(vcfOut,false);
			Assert.assertEquals(new TreeSet<>(keep), new TreeSet<>(r.getHeader().getSampleNamesInOrder()));
			r.close(); 
			}
		finally {
			support.removeTmpFiles();
			}
		}
}
