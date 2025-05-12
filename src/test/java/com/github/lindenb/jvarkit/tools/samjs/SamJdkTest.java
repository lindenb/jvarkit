package com.github.lindenb.jvarkit.tools.samjs;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;


public class SamJdkTest  {
	
	private final TestSupport support = new TestSupport();

	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		List<String> L1= support.allSamOrBams().collect(Collectors.toList());
		List<String> L2= Arrays.asList(
					"return !record.getReadUnmappedFlag() && record.getCigar().getCigarElements().stream().anyMatch(C->C.getLength()>=1000 && (C.getOperator()==CigarOperator.N || C.getOperator()==CigarOperator.D));",
					"if(record.getReadUnmappedFlag()) return true;final Cigar c=record.getCigar();if(c==null || c.numCigarElements()<2) return true; return !(c.getFirstCigarElement().getOperator().isClipping() && c.getLastCigarElement().getOperator().isClipping());"
				);
		return support.combine2(L1.stream(), L2.stream());
		}
	@Test(dataProvider="src1")
	public void test1(final String inBam,final String expr) throws IOException {
		try {
			final Path out = support.createTmpPath(".bam");
			Assert.assertEquals(new SamJdk().instanceMain(new String[] {
	        		"-o",out.toString(),
	        		"-e",expr,
	        		inBam
					}),0);
			support.assertIsValidBam(out);
			}
		finally
			{
			support.removeTmpFiles();
			}	
		}
}
