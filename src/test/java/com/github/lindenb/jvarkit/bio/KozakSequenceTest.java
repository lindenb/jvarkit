package com.github.lindenb.jvarkit.bio;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.bio.KozakSequence;


public class KozakSequenceTest {
@Test
public void test01() {
	KozakSequence k = new KozakSequence("N", 0);
	Assert.assertEquals(k.length(),KozakSequence.KOZAK_LENGTH);
	Assert.assertEquals(k.getString().length(),KozakSequence.KOZAK_LENGTH);
	Assert.assertEquals(k.getStrength(), KozakSequence.Strength.nil);
	Assert.assertEquals(k.charAt(0), 'N');
	Assert.assertEquals(k.charAt(1), 'N');
	
	k = new KozakSequence("AGGATGG",3);
	Assert.assertEquals(k.charAt(0), 'A');
	Assert.assertEquals(k.charAt(6), 'G');
	Assert.assertEquals(k.getStrength(), KozakSequence.Strength.Strong);
	k = new KozakSequence("aggatgg",3);
	Assert.assertEquals(k.getStrength(), KozakSequence.Strength.Strong);

	
	k = new KozakSequence("AGGATGC",3);
	Assert.assertEquals(k.getStrength(), KozakSequence.Strength.Moderate);

	// no atg
	k = new KozakSequence("AGGTTGG",3);
	Assert.assertEquals(k.getStrength(), KozakSequence.Strength.nil);
	}
}
