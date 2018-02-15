package com.github.lindenb.jvarkit.util.bio.fasta;

import java.io.File;
import java.io.IOException;

import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class FastaSequenceReaderTest extends TestUtils {
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new ParamCombiner().
			initList(collectAllFasta()).
			build();
	}
	
	@Test(dataProvider="src1")
	public void test01(final String inFasta) 
		throws IOException
		{
		FastaSequenceReader r=new FastaSequenceReader();
		r.getSequencesIn(new File(inFasta));
		}
	}
