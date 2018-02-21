package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class VcfToSvgTest extends TestUtils {
    @Test
    public void testVcfToSvg() throws IOException{    
    	final File svg = createTmpFile(".jeter.__SEGMENT__.svg");
    	Assert.assertTrue(svg.getName().contains("__SEGMENT__"));
    	Assert.assertEquals(0,new VcfToSvg().instanceMain(new String[]{
        		"-o",svg.getPath(),
        		"--stopAfterFirst",
        		"-k","src/test/resources/test_vcf01.knownGenes.txt.gz",
        		"src/test/resources/test_vcf01.vcf"
        		}));
    	final File svg2= new File(svg.getParentFile(),svg.getName().replaceAll("__SEGMENT__","1_861120_894679"));
    	assertIsXml(svg2);
    	Assert.assertTrue( svg2.delete());
    	Assert.assertTrue( svg.delete());
    	}

}
