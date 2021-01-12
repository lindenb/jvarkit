package com.github.lindenb.jvarkit.tools.vcfvcf;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFReader;

@AlsoTest(LauncherTest.class)
public class VcfPeekAfTest {
	private final TestSupport support =new TestSupport();

public void testACAN01()
		throws IOException
		{
		try {
			String vcfin = support.resource("test_vcf01.vcf");
			String vcfdb = vcfin;
			Path out=support.createTmpPath(".vcf");
				
			Assert.assertEquals(new VcfPeekAf().instanceMain(new String[] {
					"-o",out.toString(),
					"--peeker","ACAN",
					"--treshold","0.5",
					"--database",vcfdb,
					vcfin
					}),0);
			support.assertIsVcf(out);
			
			try(VCFReader r1= VCFReaderFactory.makeDefault().open(out,false)) {
				Assert.assertEquals(
					r1.iterator().
					stream().
					filter(V->V.getNAlleles()==2 && V.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY)).
					filter(V->V.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, 1.0) > 0.5).count(),
					0L);	
				}
			
			}
		finally {
			support.removeTmpFiles();
			}
		}

public void testGT01()
		throws IOException
		{
		try {
			String vcfin = support.resource("test_vcf01.vcf");
			String vcfdb = vcfin;
			Path out=support.createTmpPath(".vcf");
				
			Assert.assertEquals(new VcfPeekAf().instanceMain(new String[] {
					"-o",out.toString(),
					"--peeker","GT",
					"--treshold","0.5",
					"--database",vcfdb,
					vcfin
					}),0);
			support.assertIsVcf(out);
			
			try(VCFReader r1= VCFReaderFactory.makeDefault().open(out,false)) {
				Assert.assertEquals(
					r1.iterator().
					stream().
					filter(V->V.getNAlleles()==2 && V.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY)).
					filter(V->V.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, 1.0) > 0.5).count(),
					0L);	
				}
			
			}
		finally {
			support.removeTmpFiles();
			}
		}


public void testAF01()
	throws IOException
	{
	try {
		String vcfin = support.resource("test_vcf01.vcf");
		String vcfdb = vcfin;
		Path out=support.createTmpPath(".vcf");
			
		Assert.assertEquals(new VcfPeekAf().instanceMain(new String[] {
				"-o",out.toString(),
				"--peeker","AF",
				"--treshold","0.5",
				"--database",vcfdb,
				vcfin
				}),0);
		support.assertIsVcf(out);
		
		try(VCFReader r1= VCFReaderFactory.makeDefault().open(out,false)) {
			Assert.assertEquals(
				r1.iterator().
				stream().
				filter(V->V.getNAlleles()==2 && V.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY)).
				filter(V->V.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, 1.0) > 0.5).count(),
				0L);	
			}
		
		}
	finally {
		support.removeTmpFiles();
		}
	}
public void testAF02()
		throws IOException
		{
		try {
			final String fname="THEFILTER";
			String vcfin = support.resource("test_vcf01.vcf");
			String vcfdb = vcfin;
			Path out=support.createTmpPath(".vcf");
				
			Assert.assertEquals(new VcfPeekAf().instanceMain(new String[] {
					"-o",out.toString(),
					"--peeker","AF",
					"--treshold","0.5",
					"--database",vcfdb,
					"--filter",fname,
					vcfin
					}),0);
			support.assertIsVcf(out);
			
			try(VCFReader r1= VCFReaderFactory.makeDefault().open(out,false)) {
				CloseableIterator<VariantContext> iter= r1.iterator();
				while(iter.hasNext()) {
					VariantContext ctx = iter.next();
					if(ctx.getNAlleles()!=2) continue;
					if(!ctx.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY)) continue;
					double af = ctx.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, 1.0) ;
					if(af<=0.5) 
						{
						Assert.assertFalse(ctx.getFilters().contains(fname));
						}
					else
						{
						Assert.assertTrue(ctx.getFilters().contains(fname));
						}
					}
				iter.close();
				}
			
			}
		finally {
			support.removeTmpFiles();
			}
		}

}
