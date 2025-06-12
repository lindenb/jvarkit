package com.github.lindenb.jvarkit.rdf;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.rdf.ns.RDF;

public class ResourceTest {
	@Test
	public void testResource() {
		final String ns=RDF.NS;
		Assert.assertEquals(new Resource(ns,"A"), new Resource(ns,"A"));
		Assert.assertNotEquals(new Resource(ns,"A"), new Resource(ns,"a"));
		Assert.assertNotEquals(new Resource(ns,"A"), new Literal(ns+"a"));
		Assert.assertFalse(new Resource(ns,"A").isLiteral());
		Assert.assertTrue(new Resource(ns,"A").isResource());
		Assert.assertTrue(new Resource().isAnon());
		Assert.assertNotEquals(new Resource(), new Resource());
		}
}
