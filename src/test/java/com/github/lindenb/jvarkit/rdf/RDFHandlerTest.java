package com.github.lindenb.jvarkit.rdf;

import java.io.StringReader;

import com.github.lindenb.jvarkit.rdf.ns.RDF;

public class RDFHandlerTest {
	public void testRead1()  throws Exception {
		final String ns = "http://example.com/";
		final String xml="<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" +
				"<rdf:RDF xmlns:rdf=\""+RDF.NS+"\" xmlns:ex=\""+ns+"\">"+
				"<ex:Thing rdf:about=\""+ ns +"/1\">"+
				"<ex:label>LABEL</ex:label>"+
				"</ex:Thing>"+
				"</rdf:RDF>"
				;
		RDFModel m=new RDFModel();
		try(StringReader r=new StringReader(xml)) {
			RDFHandler h=new RDFHandler(STMT->m.addStatement(STMT));
			h.parse(r);
		}
	}
	

}
