package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

import htsjdk.samtools.util.IOUtil;

@AlsoTest(LauncherTest.class)
public class VcfToSvgTest {
	
	private final TestSupport support = new TestSupport();

	
    @Test
    public void testVcfToSvg() throws IOException{    
    	try {
	    	final Path svg = support.createTmpPath(".jeter.__SEGMENT__.svg");
	    	Assert.assertTrue(svg.getFileName().toString().contains("__SEGMENT__"));
	    	Assert.assertEquals(0,new VcfToSvg().instanceMain(new String[]{
	        		"-o",svg.toString(),
	        		"--stopAfterFirst",
	        		"-k",support.resource("test_vcf01.knownGenes.txt.gz"),
	        		support.resource("test_vcf01.vcf")
	        		}));
	    	final Path svg2= svg.getParent().resolve(svg.getFileName().toString().replaceAll("__SEGMENT__","1_861120_894679"));
	    	support.assertIsXml(svg2);
	    	Assert.assertTrue(Files.deleteIfExists(svg2));
	    	Assert.assertTrue(Files.deleteIfExists(svg));
	    	}
    	finally {
    		support.removeTmpFiles();
    		}
    	}
    
    @Test
    public void testVcfToSvgManifest() throws IOException{  
    	try {
	    	final Path manifest = support.createTmpPath(".MF");
	    	final String tmpTemplate = IOUtils.getDefaultTmpDir()+File.separator+"__SEGMENT__.svg";
	    	Assert.assertEquals(0,new VcfToSvg().instanceMain(new String[] {
	    			"-o",tmpTemplate,
	    			"-m",manifest.toString(),
	    			"-k",support.resource("rotavirus_rf.knowngenes.tsv.gz"),
	    			support.resource("rotavirus_rf.vcf.gz")
	    			}));
	    	Assert.assertTrue(Files.exists(manifest));
	    	Files.lines(manifest).
	    	 	map(S->S.split("[\t]")).
	    	 	map(T->T[T.length-1]).
	    	 	forEach(s->
	    		{
	    		final Path svg = Paths.get(s);
	    		support.assertIsXml(svg);
	    		IOUtil.deletePaths(svg);
	    		});
	    	}
    	finally {
    		support.removeTmpFiles();
    		}
    	}
}
