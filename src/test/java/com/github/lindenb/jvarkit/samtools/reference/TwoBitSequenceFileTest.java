package com.github.lindenb.jvarkit.samtools.reference;

import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.samtools.reference.ReferenceSequence;

public class TwoBitSequenceFileTest {
	private final TestSupport support = new TestSupport();

@Test
public void test01() throws IOException {
	String ref2bit =  support.resource("rotavirus_rf.2bit");
	Assert.assertTrue(ref2bit.endsWith(TwoBitSequenceFile.SUFFIX));
	TwoBitSequenceFile ref=null;
	try {
		ref = new TwoBitSequenceFile(ref2bit,true);
		ReferenceSequence s = ref.getSequence("xx");
		Assert.assertNull(s);
		s = ref.getSubsequenceAt("RF01",5,20);
		Assert.assertNotNull(s);
		Assert.assertNotNull(s.getBaseString());
		Assert.assertTrue(s.getBaseString().equalsIgnoreCase("attaaagctatacaAT"));
		
		s = ref.getSequence("RF11");
		Assert.assertNotNull(s);
		Assert.assertEquals(s.length(),666);
		
		}
	finally {
		if(ref!=null) ref.close();
		}
	}
}
