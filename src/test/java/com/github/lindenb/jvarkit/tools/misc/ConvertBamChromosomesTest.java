package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class ConvertBamChromosomesTest extends TestUtils {

	@Test(dataProvider="all-sam-or-bam-files")
    public void testBamRename01(final String bamin) throws IOException {
		final SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(new File(bamin));
		if(dict==null || dict.isEmpty()) return;
    	final File output = super.createTmpFile(".bam");
    	final File mapFile = super.createTmpFile("jeter.tsv");
    	final PrintWriter pw = new PrintWriter(mapFile);
    	dict.getSequences().stream().
    		map(S->S.getSequenceName()).
    	    forEach(S->{
    		pw.println(S+"\t"+(S.startsWith("chr")?S.substring(3):"chr"+S));
    		});
    	pw.flush();pw.close();
    	
        Assert.assertEquals(new ConvertBamChromosomes().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"--mapping",mapFile.getPath(),
        		bamin
        	}),0);
        assertIsValidBam(output);
    	}

}
