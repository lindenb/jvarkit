package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

@AlsoTest(LauncherTest.class)
public class ConvertVcfChromosomesTest
	{
	private final TestSupport support = new TestSupport();

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return support.toArrayArray(support.
				allVcfOrBcf().
				map(F->new Object[] {F})
				)
				;
		}
	
	@Test(dataProvider="src1")
    public void testVcfRename(final String vcfin) throws IOException {
		try {
			final SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(Paths.get(vcfin));
			if(dict==null || dict.isEmpty()) return;
	    	final Path output = support.createTmpPath(".vcf");
	    	final Path mapFile = support.createTmpPath("jeter.tsv");
	    	final PrintWriter pw = new PrintWriter(Files.newBufferedWriter(mapFile));
	    	dict.getSequences().stream().
	    		map(S->S.getSequenceName()).
	    	    forEach(S->{
	    		pw.println(S+"\t"+(S.startsWith("chr")?S.substring(3):"chr"+S));
	    		});
	    	pw.flush();pw.close();
	    	
	        Assert.assertEquals(0,new ConvertVcfChromosomes().instanceMain(new String[]{
	        		"-o",output.toString(),
	        		"--mapping",mapFile.toString(),
	        		vcfin
	        	}));
	        support.assertIsVcf(output);
	        Assert.assertTrue( support.variantStream(output).noneMatch(V->dict.getSequence(V.getContig())!=null));
			}
		finally {
				support.removeTmpFiles();
			}
		}
}
