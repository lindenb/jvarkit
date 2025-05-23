package com.github.lindenb.jvarkit.tools.vcf2intervals;

import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class VcfToIntervalsTest {
	@Test
	public void testNVariants() throws Exception {
		final TestSupport support=new TestSupport();
		try {
			Path bed= support.createTmpPath(".bed");
			Assert.assertEquals(new VcfToIntervals().instanceMain(new String[] {
					"-N","10","--bed",
					"-o",bed.toString(),
					support.resource("rotavirus_rf.vcf.gz")})
					,0);
			
			support.assertIsBed(bed);
			}
		catch(Throwable err) {
			Assert.fail("testNVariants", err);
			}
		finally {
			support.removeTmpFiles();
			}
		}
	
	@Test
	public void testDistance() throws Exception {
		final TestSupport support=new TestSupport();
		try {
			Path bed= support.createTmpPath(".bed");
			Assert.assertEquals(new VcfToIntervals().instanceMain(new String[] {
					"--distance","500",
					"--bed",
					"-o",bed.toString(),
					support.resource("rotavirus_rf.vcf.gz")})
					,0);
			
			support.assertIsBed(bed);
			}
		catch(Throwable err) {
			Assert.fail("testDistance", err);
			}
		finally {
			support.removeTmpFiles();
			}
		}
}
