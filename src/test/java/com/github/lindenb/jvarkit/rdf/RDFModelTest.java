package com.github.lindenb.jvarkit.rdf;

import java.io.IOException;
import java.nio.file.Path;

import javax.xml.stream.XMLStreamException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class RDFModelTest {
@Test
public void testAdd() throws IOException,XMLStreamException {
	final RDFModel m=new RDFModel();
	final Resource r1 = new Resource("http://example.com/","X");
	for(int i=0;i< 10000;i++) {
		m.addStatement(r1,r1, "Hello");
		m.addStatement(r1,r1, "World");
		Assert.assertEquals(m.size(), 2);
		Assert.assertEquals(m.findMatching(null, null, null).count(),2L);
		Assert.assertEquals(m.findMatching(r1, null, null).count(),2L);
		Assert.assertEquals(m.findMatching(null,r1,null).count(),2L);
		Assert.assertEquals(m.findMatching(null,null,new Literal("Hello")).count(),1L);
		Assert.assertEquals(m.findMatching(r1,null,new Literal("Hello")).count(),1L);
		Assert.assertEquals(m.findMatching(r1,r1,new Literal("Hello")).count(),1L);
		Assert.assertEquals(m.findMatching(r1,r1,new Literal("None")).count(),0L);
		
		
		}
	final TestSupport support=new TestSupport();
	try {
		Path p = support.createTmpPath(".rdf");
		m.writeXml(p);
		}
	finally {
		support.removeTmpFiles();
		}
	m.removeMatching(null,null,new Literal("Hello"));
	Assert.assertEquals(m.size(),1);
	m.removeMatching(r1,null,null);
	Assert.assertEquals(m.size(),0);
	}

}
