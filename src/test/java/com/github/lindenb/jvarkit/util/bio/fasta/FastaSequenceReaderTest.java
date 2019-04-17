package com.github.lindenb.jvarkit.util.bio.fasta;

import java.io.File;
import java.io.IOException;

import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class FastaSequenceReaderTest {
	
	private final TestSupport support =new TestSupport();

	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return support.toArrayArray(support.allFasta().map(F->new Object[] {F}));
	}
	
	@Test(dataProvider="src1")
	public void test01(final String inFasta) 
		throws IOException
		{
		try {
			final FastaSequenceReader r=new FastaSequenceReader();
			r.getSequencesIn(new File(inFasta));
		} finally {
			support.removeTmpFiles();
			}
		}
	}
