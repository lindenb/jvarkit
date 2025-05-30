package com.github.lindenb.jvarkit.rdf;

import org.testng.Assert;
import org.testng.annotations.Test;

public class LiteralTest {
@Test
public void testLiteral() {
	Assert.assertEquals(new Literal("A"), new Literal("A"));
	Assert.assertNotEquals(new Literal("A"), new Literal("a"));
	Assert.assertTrue(new Literal("A").isLiteral());
	Assert.assertFalse(new Literal("A").isResource());
	}
}
