package com.github.lindenb.jvarkit.bio;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.bio.Rebase;

public class RebaseTest {
@Test
public void test() {
	final Rebase rebase = Rebase.createDefaultRebase();
	Assert.assertNotEquals(0,rebase.size());
	Rebase.Enzyme enz = rebase.getEnzymeByName("hello");
	Assert.assertNull(enz);
	enz = rebase.getEnzymeByName("EcoRI");
	Assert.assertNotNull(enz);
	Assert.assertEquals(enz.getName(), "EcoRI");
	Assert.assertEquals(enz.getDecl(), "G^AATTC");
	Assert.assertEquals(enz.getWeight(), 6f);
	Assert.assertEquals(enz.getBases(), "GAATTC");
	Assert.assertTrue(enz.isPalindromic());
	}
}
