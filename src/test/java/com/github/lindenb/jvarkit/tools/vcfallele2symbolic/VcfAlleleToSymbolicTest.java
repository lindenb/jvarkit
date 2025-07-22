package com.github.lindenb.jvarkit.tools.vcfallele2symbolic;

import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class VcfAlleleToSymbolicTest {
	@Test
	public void testSymbolic1() throws Exception {
		TestSupport support=new TestSupport();
		Path vin = null;
        Path vout = null;
        try {
            final String vcf=     
            "##fileformat=VCFv4.2\n"+
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"+
            "1\t1\t.\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tC\t0\tPASS\t.\n" + 
            "1\t2\t.\tA\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA,CCCCCCCCCCCCCCCCCCCC\t0\tPASS\t.\n"  
            ;
            vin= Files.createTempFile("tmp", ".vcf");
            vout= Files.createTempFile("tmp", ".vcf");
            Files.writeString(vin, vcf);
            
            Assert.assertEquals(0,new VcfAlleleToSymbolic().instanceMain(new String[] {"-n","10","--keep-ref","-o",vout.toString(),vin.toString()}));
			}
		finally {
			support.removeTmpFiles();
		}
	}

}
