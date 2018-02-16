package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class ConvertVcfChromosomesTest extends TestUtils 
	{
	@Test(dataProvider="all-vcf-files")
    public void testVcfRename(final String vcfin) throws IOException {
		final SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(new File(vcfin));
		if(dict==null || dict.isEmpty()) return;
    	final File output = super.createTmpFile(".vcf");
    	final File mapFile = super.createTmpFile("jeter.tsv");
    	final PrintWriter pw = new PrintWriter(mapFile);
    	dict.getSequences().stream().
    		map(S->S.getSequenceName()).
    	    forEach(S->{
    		pw.println(S+"\t"+(S.startsWith("chr")?S.substring(3):"chr"+S));
    		});
    	pw.flush();pw.close();
    	
        Assert.assertEquals(0,new ConvertVcfChromosomes().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"--mapping",mapFile.getPath(),
        		vcfin
        	}));
        assertIsVcf(output);
        Assert.assertTrue(variantStream(output).noneMatch(V->dict.getSequence(V.getContig())!=null));
    	}
}
