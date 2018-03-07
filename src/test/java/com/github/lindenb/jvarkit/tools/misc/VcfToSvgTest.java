package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.samtools.util.IOUtil;

public class VcfToSvgTest extends TestUtils {
    @Test
    public void testVcfToSvg() throws IOException{    
    	final File svg = createTmpFile(".jeter.__SEGMENT__.svg");
    	Assert.assertTrue(svg.getName().contains("__SEGMENT__"));
    	Assert.assertEquals(0,new VcfToSvg().instanceMain(new String[]{
        		"-o",svg.getPath(),
        		"--stopAfterFirst",
        		"-k",SRC_TEST_RESOURCE + "/test_vcf01.knownGenes.txt.gz",
        		SRC_TEST_RESOURCE + "/test_vcf01.vcf"
        		}));
    	final File svg2= new File(svg.getParentFile(),svg.getName().replaceAll("__SEGMENT__","1_861120_894679"));
    	assertIsXml(svg2);
    	Assert.assertTrue( svg2.delete());
    	Assert.assertTrue( svg.delete());
    	}
    @Test
    public void testVcfToSvgManifest() throws IOException{    
    	final File manifest = createTmpFile(".MF");
    	final String tmpTemplate = IOUtils.getDefaultTmpDir()+File.separator+"__SEGMENT__.svg";
    	Assert.assertEquals(0,new VcfToSvg().instanceMain(new String[] {
    			"-o",tmpTemplate,
    			"-m",manifest.getPath(),
    			"-k",SRC_TEST_RESOURCE + "/rotavirus_rf.knowngenes.tsv.gz",
        		SRC_TEST_RESOURCE + "/rotavirus_rf.vcf.gz"
    			}));
    	Assert.assertTrue(manifest.exists());
    	 IOUtil.slurpLines(manifest).stream().
    	 	map(S->S.split("[\t]")).
    	 	map(T->T[T.length-1]).
    	 	forEach(s->
    		{
    		final File svg = new File(s);
    		assertIsXml(svg);
    		Assert.assertTrue(svg.delete());
    		});
    	}
}
