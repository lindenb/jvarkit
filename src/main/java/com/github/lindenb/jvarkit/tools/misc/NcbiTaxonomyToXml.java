/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.samtools.util.CloserUtil;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.command.Command;

public class NcbiTaxonomyToXml extends AbstractNcbiTaxonomyToXml
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(NcbiTaxonomyToXml.class);

	
	 @Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	
		
	private static class Node
		{
		int id=0;
		int parent_id=-1;
		Set<Node> children=new HashSet<Node>();;
		String common_name=null;
		String scientific_name=null;
		String rank=null;
		
		private Node(int id)
			{
			this.id=id;
			common_name=String.valueOf(this.id);
			scientific_name=common_name;
			}	
		@Override
		public int hashCode()
			{
			return id;
			}
		@Override
		public boolean equals(Object obj)
			{
			if(obj==this) return true;
			return Node.class.cast(obj).id==this.id;
			}
		
		
		void writeXml(XMLStreamWriter w) throws XMLStreamException
			{
			w.writeStartElement("Taxon");
			w.writeAttribute("id", String.valueOf(this.id));
			w.writeStartElement("TaxId");
			w.writeCharacters(String.valueOf(this.id));
			w.writeEndElement();
			if(this.rank!=null && !this.rank.equals("no rank"))
				{
				w.writeStartElement("Rank");
				w.writeCharacters(this.rank);
				w.writeEndElement();
				}
			w.writeStartElement("ScientificName");
			w.writeCharacters(this.scientific_name);
			w.writeEndElement();
			if(parent_id!=-1)
				{
				w.writeStartElement("ParentTaxId");
				w.writeCharacters(String.valueOf(parent_id));
				w.writeEndElement();
				}

			w.writeStartElement("TaxaSet");
			for(Node c:this.children)
				{
				c.writeXml(w);
				}
			w.writeEndElement();
			w.writeEndElement();
			}
		
		}
	public  class MyCommand extends AbstractNcbiTaxonomyToXml.AbstractNcbiTaxonomyToXmlCommand
	 	{
		private Map<Integer,Node> id2node=new HashMap<Integer,Node>();
	
		
		private Node getNodeById(int id)
			{
			Node node=id2node.get(id);
			if(node==null)
				{
				node=new Node(id);
				id2node.put(id,node);
				}
			return node;
			}
		
		private void readNames(BufferedReader in) throws IOException
			{
			Pattern delim=Pattern.compile("\t\\|(\t)?");
			String line;
			while((line=in.readLine())!=null)
				{
				String tokens[]=delim.split(line);
				int id= Integer.parseInt(tokens[0].trim());
				Node node=this.id2node.get(id);
				if(node==null) continue;
				
				if(tokens[3].equals("common name"))
					{
					node.common_name=tokens[1].trim();
					}
				else if(tokens[3].equals("scientific name"))
					{
					node.scientific_name=tokens[1].trim();
					}
				}
			}	
		private void readNodes(BufferedReader in) throws IOException
			{
			Pattern delim=Pattern.compile("\t\\|(\t)?");
			String line;
			while((line=in.readLine())!=null)
				{
				String tokens[]=delim.split(line);
				
				Node node = getNodeById(Integer.parseInt( tokens[0].trim()));
				Node parent= getNodeById(Integer.parseInt( tokens[1].trim()));
				if(parent!=node)
					{
					parent.children.add(node);
					node.parent_id=parent.id;
					}
				node.rank=tokens[2].trim();
				}
			}	
		
		@Override
		protected Collection<Throwable> call(File baseDir) throws Exception
			{
			PrintStream out=null;
			try
				{
				File inputFile=new File(baseDir,"nodes.dmp");
				LOG.info("Reading Nodes "+inputFile);
				BufferedReader in=IOUtils.openFileForBufferedReading(inputFile);
				readNodes(in);
				in.close();
				
				inputFile=new File(baseDir,"names.dmp");
				LOG.info("Reading Names "+inputFile);
				in=IOUtils.openFileForBufferedReading(inputFile);
				readNames(in);
				in.close();
				
				Node root= this.id2node.get(root_id);
				if(root==null)
					{
					return wrapException("Cannot get node id."+root_id);
					}
				
				out = openFileOrStdoutAsPrintStream();
				XMLOutputFactory xof=XMLOutputFactory.newFactory();
				XMLStreamWriter w=xof.createXMLStreamWriter(out, "UTF-8");
				w.writeStartDocument("UTF-8", "1.0");
				w.writeStartElement("TaxaSet");
				w.writeComment(getName()+" by "+getAuthorName());
				//w.writeComment(getP());
				root.writeXml(w);
				w.writeEndElement();//TaxaSet
		
				w.writeEndDocument();
				w.flush();
				out.flush();
				return RETURN_OK;
				}
			catch(Exception err)
				{
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(out);
				}
			}
	 	}
	
	public static void main(String[] args)
		{
		new NcbiTaxonomyToXml().instanceMainWithExit(args);
		}
	}
