package com.github.lindenb.jvarkit.tools.misc;


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
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;

@Program(name="ncbitaxonomy2xml",
	description="Dump NCBI taxonomy tree as a hierarchical XML document",
	keywords={"taxonomy","ncbi","xml"}
	)
public class NcbiTaxonomyToXml extends Launcher
	{
	private static final Logger LOG = Logger.build(NcbiTaxonomyToXml.class).make();
	
	@Parameter(names={"-o","--out"},description="Output file or stdout")
	private File outputFile = null;
	@Parameter(names={"-r","--root"},description="NCBI taxon root id.")
	private int root_id=1;

	
	private final Map<Integer,Node> id2node=new HashMap<Integer,Node>();

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
		
		
		void writeXml(final XMLStreamWriter w) throws XMLStreamException
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
	
	private NcbiTaxonomyToXml()
		{
		
		}
	
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
	public int doWork(List<String> args) {
		OutputStream out=null;
		File baseDir=null;
		try
			{
			final String baseDirStr = super.oneFileOrNull(args);
			if(baseDirStr==null)
				{
				baseDir=null;
				}
			else 
				{
				baseDir=new File(baseDirStr);
				}
			
			
			
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
				LOG.error("Cannot get node id."+root_id);
				return -1;
				}
			
			out= super.openFileOrStdoutAsStream(outputFile);
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			XMLStreamWriter w=xof.createXMLStreamWriter(out, "UTF-8");
			w.writeStartDocument("UTF-8", "1.0");
			w.writeStartElement("TaxaSet");
			w.writeComment(getProgramCommandLine());
			root.writeXml(w);
			w.writeEndElement();//TaxaSet
	
			w.writeEndDocument();
			w.flush();
			out.flush();
			out.close();out=null;
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
			}
		}
	public static void main(String[] args)
		{
		new NcbiTaxonomyToXml().instanceMainWithExit(args);
		}
	}
