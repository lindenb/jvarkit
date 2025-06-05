package com.github.lindenb.jvarkit.rdf;


import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.rdf.ns.RDF;

public class StatementTest {
	@Test
	public void testSmts() {
		final String ns=RDF.NS;
		final Resource r1 = new Resource(ns,"A");
		final Resource r2 = new Resource(ns,"B");
		final Resource p = new Resource(ns,"C");
		final RDFNode o1 = new Literal("Hello");
		
		final Statement stmt1 = new Statement(r1, p, o1);
		final Statement stmt2 = new Statement(r1, p, o1);
		final Statement stmt3 = new Statement(r2, p, o1);
		Assert.assertEquals(stmt1,stmt2);
		Assert.assertNotEquals(stmt1,stmt3);
		Assert.assertEquals(stmt1.getSubject(),stmt2.getSubject());
		Assert.assertEquals(stmt1.getPredicate(),stmt2.getPredicate());
		Assert.assertEquals(stmt1.getObject(),stmt2.getObject());
		}
}
