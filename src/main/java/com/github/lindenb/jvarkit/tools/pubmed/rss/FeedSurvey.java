package com.github.lindenb.jvarkit.tools.pubmed.rss;

import java.io.StringReader;
import java.util.List;

import javax.xml.parsers.DocumentBuilder;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.xml.sax.InputSource;

public class FeedSurvey {
	DocumentBuilder db;

	private class Entity {
		
		}
	private class GeneEntity extends Entity{
		public boolean find(List<String> words) {
			return false;
			}
		}
	private class DiseaseEntity extends Entity{
		
		}
	
	private void findEntities(final String s,final List<Entity> entities) {
		
		}
	
	private void item(Element item) throws Exception {
		for(Node c1=item.getFirstChild();c1!=null;c1=c1.getNextSibling()) {
			Element content_encodded = Element.class.cast(c1);
			String c= content_encodded.getTextContent();
			try(StringReader sr=new StringReader(c)) {
				Document textDom = db.parse(new InputSource(sr));
				Element div = textDom.getDocumentElement();
				}
			}
		}
	
	private void load(String feedURL) throws Exception {
		Document dom=db.parse(feedURL);
		Element rss = dom.getDocumentElement();
		for(Node c1=rss.getFirstChild();c1!=null;c1=c1.getNextSibling()) {
			Element channel = Element.class.cast(c1);
			for(Node c2=c1.getFirstChild();c2!=null;c2=c2.getNextSibling()) {
				Element item = Element.class.cast(c2);
				
				}
			}
		}
}
