package com.github.lindenb.jvarkit.tools.taxonomy;


/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
import java.io.PrintWriter;
import java.nio.file.Path;
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
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.log.Logger;

import htsjdk.samtools.util.CloserUtil;

/*
 * 
BEGIN_DOC
## Example

producing taxonomy.xml using make:

```bash
taxonomy.xml: nodes.dmp names.dmp
	java -jar dist/ncbitaxonomy2xml.jar . | xmllint --format - > $@

nodes.dmp : taxdump.tar.gz
	tar xvfz $< $@
names.dmp :taxdump.tar.gz
	tar xvfz $< $@

taxdump.tar.gz:
	curl -o $@  "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"

```

output:
```xml
<?xml version="1.0" encoding="UTF-8"?>
<TaxaSet>
  <!--NcbiTaxonomyToXml by Pierre Lindenbaum PhD.-->
  <!--.-->
  <Taxon id="1">
    <TaxId>1</TaxId>
    <ScientificName>root</ScientificName>
    <TaxaSet>
      <Taxon id="131567">
        <TaxId>131567</TaxId>
        <ScientificName>cellular organisms</ScientificName>
        <ParentTaxId>1</ParentTaxId>
        <TaxaSet>
          <Taxon id="2">
            <TaxId>2</TaxId>
            <Rank>superkingdom</Rank>
            <ScientificName>Bacteria</ScientificName>
            <ParentTaxId>131567</ParentTaxId>
            <TaxaSet>
              <Taxon id="51290">
                <TaxId>51290</TaxId>
                <Rank>superphylum</Rank>
                <ScientificName>Chlamydiae/Verrucomicrobia group</ScientificName>
                <ParentTaxId>2</ParentTaxId>
                <TaxaSet>
                  <Taxon id="256845">
                    <TaxId>256845</TaxId>
                    <Rank>phylum</Rank>
                    <ScientificName>Lentisphaerae</ScientificName>
                    <ParentTaxId>51290</ParentTaxId>
                    <TaxaSet>
                      <Taxon id="278094">
                        <TaxId>278094</TaxId>
                        <ScientificName>environmental samples</ScientificName>
                        <ParentTaxId>256845</ParentTaxId>
                        <TaxaSet>
                          <Taxon id="278095">
                            <TaxId>278095</TaxId>
                            <Rank>species</Rank>
                            <ScientificName>uncultured Lentisphaerae bacterium</ScientificName>
                            <ParentTaxId>278094</ParentTaxId>
                            <TaxaSet/>
                          </Taxon>
(...)
```

END_DOC
*/

@Program(name="ncbitaxonomy2xml",
	description="Dump NCBI taxonomy tree as a hierarchical XML document or as a table",
	keywords={"taxonomy","ncbi","xml"},
	creationDate = "20120320",
	modificationDate = "20240320",
	biostars=10327,
	jvarkit_amalgamion = true,
	menu="Utilities"
	)
public class NcbiTaxonomyToXml extends Launcher
	{
	private static final Logger LOG = Logger.of(NcbiTaxonomyToXml.class);
	
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-r","--root"},description="NCBI taxon root id.")
	private int root_id=1;
	@Parameter(names={"-t","--tabular"},description="produce tabular output instead of xml")
	private boolean tabular_output=false;

	
	private final Map<Integer,Node> id2node=new HashMap<Integer,Node>();

	private static class Node
		{
		int id=0;
		int parent_id=-1;
		final Set<Node> children=new HashSet<Node>();;
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
			return Integer.hashCode(id);
			}
		@Override
		public boolean equals(Object obj)
			{
			if(obj==this) return true;
			return Node.class.cast(obj).id==this.id;
			}
		
		void writeTabular(final PrintWriter w) {
			w.print(String.valueOf(this.id));
			w.print("\t");
			w.print(this.parent_id!=-1?String.valueOf(this.parent_id):".");
			w.print("\t");
			w.print(this.scientific_name);
			w.print("\t");
			w.print(this.rank!=null && !this.rank.equals("no rank")?this.rank:".");
			w.println();
			
			for(Node c:this.children)
				{
				c.writeTabular(w);
				}
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
	
	public NcbiTaxonomyToXml()
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
		final Pattern delim=Pattern.compile("\t\\|(\t)?");
		String line;
		while((line=in.readLine())!=null)
			{
			final String tokens[]=delim.split(line);
			int id= Integer.parseInt(tokens[0].trim());
			final Node node=this.id2node.get(id);
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
		final Pattern delim=Pattern.compile("\t\\|(\t)?");
		String line;
		while((line=in.readLine())!=null)
			{
			final String tokens[]=delim.split(line);
			
			final Node node = getNodeById(Integer.parseInt( tokens[0].trim()));
			final Node parent= getNodeById(Integer.parseInt( tokens[1].trim()));
			if(parent!=node)
				{
				parent.children.add(node);
				node.parent_id=parent.id;
				}
			node.rank=tokens[2].trim();
			}
		}	
	
	
	@Override
	public int doWork(final List<String> args) {
		try
			{
			final File baseDir;
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
			try(BufferedReader in=IOUtils.openFileForBufferedReading(inputFile)) {
				readNodes(in);
				}
			
			inputFile=new File(baseDir,"names.dmp");
			LOG.info("Reading Names "+inputFile);
			try(BufferedReader in=IOUtils.openFileForBufferedReading(inputFile)) {
				readNames(in);
				}
			
			final Node root= this.id2node.get(root_id);
			if(root==null)
				{
				LOG.error("Cannot get node id."+root_id);
				return -1;
				}
			if(tabular_output) {
				try(PrintWriter out= super.openPathOrStdoutAsPrintWriter(outputFile)) {
					out.println("#taxon\tparent\tname\trank");
					root.writeTabular(out);
					out.flush();
					}
				}
			else
				{
				try(OutputStream out= super.openPathOrStdoutAsStream(outputFile)) {
					final XMLOutputFactory xof=XMLOutputFactory.newFactory();
					final XMLStreamWriter w=xof.createXMLStreamWriter(out, "UTF-8");
					w.writeStartDocument("UTF-8", "1.0");
					w.writeStartElement("TaxaSet");
					w.writeComment(getProgramCommandLine());
					root.writeXml(w);
					w.writeEndElement();//TaxaSet
			
					w.writeEndDocument();
					w.flush();
					out.flush();
					}
				}
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		}
	public static void main(final String[] args)
		{
		new NcbiTaxonomyToXml().instanceMainWithExit(args);
		}
	}
