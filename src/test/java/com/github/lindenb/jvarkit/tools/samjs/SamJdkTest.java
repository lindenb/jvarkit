package com.github.lindenb.jvarkit.tools.samjs;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;


public class SamJdkTest extends TestUtils {
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new ParamCombiner().
			initList(collectAllSamOrBam()).
			product(
				"return !record.getReadUnmappedFlag() && record.getCigar().getCigarElements().stream().anyMatch(C->C.getLength()>=1000 && (C.getOperator()==CigarOperator.N || C.getOperator()==CigarOperator.D));",
				"if(record.getReadUnmappedFlag()) return true;final Cigar c=record.getCigar();if(c==null || c.numCigarElements()<2) return true; return !(c.getFirstCigarElement().getOperator().isClipping() && c.getLastCigarElement().getOperator().isClipping());").
			build();
		}
	@Test(dataProvider="src1")
	public void test1(final String inBam,final String expr) throws IOException {
		final File out = createTmpFile(".bam");
		Assert.assertEquals(0,new SamJdk().instanceMain(newCmd().add(
        		"-o",out.getPath(),
        		"-e",expr,
        		inBam
        		).
				make()));
		assertIsValidBam(out);
		}
}
