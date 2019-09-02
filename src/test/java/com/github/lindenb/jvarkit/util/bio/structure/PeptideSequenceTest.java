package com.github.lindenb.jvarkit.util.bio.structure;

import org.testng.Assert;
import org.testng.annotations.Test;

public class PeptideSequenceTest {
@Test
public void test01() {
	String rna="Atg";
	PeptideSequence<String> pep = PeptideSequence.of(rna);
	Assert.assertNotNull(pep);
	Assert.assertEquals(1, pep.length());
	Assert.assertEquals(pep.charAt(0),'M');
	Assert.assertEquals(pep.toString(),"M");
	}
}
