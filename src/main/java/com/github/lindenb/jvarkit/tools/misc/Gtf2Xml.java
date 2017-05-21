/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.samtools.util.CloserUtil;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StreamTokenizer;
import java.io.StringReader;
import java.net.URLDecoder;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
/*
BEGIN_DOC


## Example

```bash
$ curl  "ftp://ftp.ensembl.org/pub/release-81/gff3/homo_sapiens/Homo_sapiens.GRCh38.81.gff3.gz" | gunzip -c |\
 head -n 10 | java -jar dist/gtf2xml.jar | xmllint --format -

```
output: 
```xml
<?xml version="1.0" encoding="UTF-8"?>
<gtf gff-version="3" genome-build="GRCh38.p3" genome-version="GRCh38" genome-date="2013-12" genome-build-accession="NCBI:GCA_000001405.18" genebuild-last-updated="2015-06">
  <feature chrom="1" start="11869" end="14409" strand="+" source="havana" type="gene">
    <attributes>
      <ID>gene:ENSG00000223972</ID>
      <Name>DDX11L1</Name>
      <biotype>transcribed_unprocessed_pseudogene</biotype>
      <description>DEAD/H (Asp-Glu-Ala-Asp/His) box helicase 11 like 1 [Source:HGNC Symbol;Acc:HGNC:37102]</description>
      <gene_id>ENSG00000223972</gene_id>
      <havana_gene>OTTHUMG00000000961</havana_gene>
      <havana_version>2</havana_version>
      <logic_name>havana</logic_name>
      <version>5</version>
    </attributes>
  </feature>
  <feature chrom="1" start="14404" end="29570" strand="-" source="havana" type="gene">
    <attributes>
      <ID>gene:ENSG00000227232</ID>
      <Name>WASH7P</Name>
      <biotype>unprocessed_pseudogene</biotype>
      <description>WAS protein family homolog 7 pseudogene [Source:HGNC Symbol;Acc:HGNC:38034]</description>
      <gene_id>ENSG00000227232</gene_id>
      <havana_gene>OTTHUMG00000000958</havana_gene>
      <havana_version>1</havana_version>
      <logic_name>havana</logic_name>
      <version>5</version>
    </attributes>
  </feature>
  <feature chrom="1" start="17369" end="17436" strand="-" source="ensembl" type="miRNA_gene">
    <attributes>
      <ID>gene:ENSG00000278267</ID>
      <Name>MIR6859-1</Name>
      <biotype>miRNA</biotype>
      <description>microRNA 6859-1 [Source:HGNC Symbol;Acc:HGNC:50039]</description>
      <gene_id>ENSG00000278267</gene_id>
      <logic_name>ncrna</logic_name>
      <version>1</version>
    </attributes>
  </feature>
  <feature chrom="1" start="29554" end="31109" strand="+" source="havana" type="lincRNA_gene">
    <attributes>
      <ID>gene:ENSG00000243485</ID>
      <Name>RP11-34P13.3</Name>
      <biotype>lincRNA</biotype>
      <gene_id>ENSG00000243485</gene_id>
      <havana_gene>OTTHUMG00000000959</havana_gene>
      <havana_version>2</havana_version>
      <logic_name>havana</logic_name>
      <version>3</version>
    </attributes>
  </feature>
  <attributes>
    <attribute>havana_gene</attribute>
    <attribute>Name</attribute>
    <attribute>havana_version</attribute>
    <attribute>logic_name</attribute>
    <attribute>description</attribute>
    <attribute>biotype</attribute>
    <attribute>ID</attribute>
    <attribute>gene_id</attribute>
    <attribute>version</attribute>
  </attributes>
  <types>
    <type>miRNA_gene</type>
    <type>gene</type>
    <type>lincRNA_gene</type>
  </types>
  <sources>
    <source>ensembl</source>
    <source>havana</source>
  </sources>
  <dict>
    <chrom name="1" length="31109"/>
  </dict>
</gtf>
```

 
END_DOC
 */
@Program(name="gtf2xml",
	description="Convert GTF to XML",
	keywords={"xml","gtf"})
public class Gtf2Xml extends Launcher{
	private static final Logger LOG = Logger.build(FixVCF.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-T","--tmpDir"},description="mp directory")
	private File tmpDir = new File(System.getProperty("java.io.tmpdir"));

	private abstract class GffCodec
		{
		Map<String,Long> seqdict=new LinkedHashMap<>();
		Set<String> att_keys=new HashSet<>();
		Set<String> sources=new HashSet<>();
		Set<String> types=new HashSet<>();
		protected Pattern tab=Pattern.compile("[\t]");
		protected Pattern semicolon=Pattern.compile("[;]");
		void write(XMLStreamWriter w,String line) throws XMLStreamException,IOException
			{
			String tokens[] = this.tab.split(line);
			if(tokens.length<8) throw new IOException("Expected at least 8 columns in "+line+" got "+tokens.length);
			w.writeStartElement("feature");
			
			w.writeAttribute("chrom", tokens[0]);
			w.writeAttribute("start", tokens[3]);
			w.writeAttribute("end", tokens[4]);
			
			Long contifLength  = seqdict.get(tokens[0]);
			if(contifLength==null) contifLength=0L;
			seqdict.put(tokens[0],Math.max(Long.parseLong(tokens[4]), contifLength));
			
			if(!tokens[5].equals("."))  w.writeAttribute("score", tokens[5]);
			if(!tokens[6].equals("."))  w.writeAttribute("strand", tokens[6]);
			if(!tokens[1].equals("."))
				{
				w.writeAttribute("source", tokens[1]);
				sources.add(tokens[1]);
				}
			if(!tokens[2].equals("."))
				{
				w.writeAttribute("type", tokens[2]);
				types.add(tokens[2]);
				}
			if(!tokens[7].equals(".")) w.writeAttribute("phase", tokens[7]);
			if(!tokens[8].equals("."))
				{
				w.writeStartElement("attributes");
				writeAttributes(w,tokens[8]);
				w.writeEndElement();
				}
			w.writeEndElement();
			w.writeCharacters("\n");
			}
		abstract void writeAttributes(XMLStreamWriter w,String attrString)  throws XMLStreamException,IOException;
			
		
		}
	
	private class DefaultCodec extends GffCodec
		{
		@Override
		void writeAttributes(XMLStreamWriter w,String attrString) throws XMLStreamException,IOException
			{
			StreamTokenizer st=new StreamTokenizer(new StringReader(attrString));
			st.wordChars('_', '_');
			String key=null;
			while(st.nextToken() != StreamTokenizer.TT_EOF)
				{
				String s=null;
				switch(st.ttype)
					{
					case StreamTokenizer.TT_NUMBER: s=String.valueOf(st.nval);break;
					case '"': case '\'' : case StreamTokenizer.TT_WORD: s=st.sval;break;
					case ';':break;
					default:break;
					}
				if(s==null) continue;
				if(key==null)
					{
					key=s;
					this.att_keys.add(key);
					}
				else 
					{
					w.writeStartElement(key);
					w.writeCharacters(s);
					w.writeEndElement();
					key=null;
					}
				}
			
			}
		}
	
	private class Gtf3Codec extends GffCodec
		{
		private String unescape(String s) throws IOException
			{
			return URLDecoder.decode(s, "UTF-8");
			}
		@Override
		void writeAttributes(XMLStreamWriter w,String attrString) throws XMLStreamException,IOException
			{
			String atts[]  = semicolon.split(attrString);
			for(String att: atts)
				{
				if(att.isEmpty()) continue;
				int eq=att.indexOf("=");
				if(eq<=0) throw new IOException("bad att "+att+" in "+attrString);
				String key = att.substring(0,eq);
				String value = att.substring(eq+1);
				w.writeStartElement(key);
				w.writeCharacters(unescape(value) );
				w.writeEndElement();
				this.att_keys.add(key);
				}
			}
		}
	public Gtf2Xml() {
		}
	
	@Override
	public int doWork(List<String> args) {

		BufferedReader r=null;
		XMLStreamWriter w=null;
		FileWriter fw=null;
		try {
			String inputName=oneFileOrNull(args);
			r = super.openBufferedReader(inputName);
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			if(this.outputFile==null)
				{
				w = xof.createXMLStreamWriter(stdout(), "UTF-8");
				}
			else
				{
				w = xof.createXMLStreamWriter((fw=new FileWriter(this.outputFile)));
				}
			GffCodec codec = new DefaultCodec();
			String headerLine=null;
			w.writeStartDocument("UTF-8","1.0");
			w.writeStartElement("gtf");
			while((headerLine=r.readLine())!=null)
				{
				if(!headerLine.startsWith("#")) break;
				int ws = headerLine.indexOf(' ');
				if(ws==-1) continue; //??
				
				if(headerLine.startsWith("##gff-version "))
					{
					String version=headerLine.substring(ws+1).trim();
					LOG.info("version "+version);
					if(version.equals("3"))
						{
						codec = new Gtf3Codec();
						}
					w.writeAttribute("gff-version", version);
					}
				else if(headerLine.startsWith("#!"))
					{
					w.writeAttribute(
							headerLine.substring(2, ws),
							headerLine.substring(ws+1).trim());
					}
				else
					{
					LOG.warning("ignoring "+headerLine);
					}
				}
			w.writeCharacters("\n");
			
			
		
			String line=null;
			for(;;)
				{
				line=(headerLine!=null?headerLine:r.readLine());
				headerLine=null;
				if(line==null) break;
				if(line.isEmpty()) continue;
				codec.write(w,line);
				}
			
			
			w.writeStartElement("attributes");
			for(String k : codec.att_keys)
				{
				w.writeStartElement("attribute");
				w.writeCharacters(k);
				w.writeEndElement();
				w.writeCharacters("\n");
				}
			w.writeEndElement();
			w.writeCharacters("\n");
			
			w.writeStartElement("types");
			for(String k : codec.types)
				{
				w.writeStartElement("type");
				w.writeCharacters(k);
				w.writeEndElement();
				w.writeCharacters("\n");
				}
			w.writeEndElement();
			w.writeCharacters("\n");
			
			
			w.writeStartElement("sources");
			for(String k : codec.sources)
				{
				w.writeStartElement("source");
				w.writeCharacters(k);
				w.writeEndElement();
				w.writeCharacters("\n");
				}
			w.writeEndElement();
			w.writeCharacters("\n");

			w.writeStartElement("dict");
			for(String k : codec.seqdict.keySet())
				{
				w.writeEmptyElement("chrom");
				w.writeAttribute("name",k);
				w.writeAttribute("length", String.valueOf(codec.seqdict.get(k)));
				w.writeCharacters("\n");
				}
			w.writeEndElement();
			w.writeCharacters("\n");

			
			w.writeEndElement();
			w.writeEndDocument();
			w.flush();
			return 0;
			}
		catch (Exception e) {
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);
			CloserUtil.close(fw);
			}
		}

	
	public static void main(String[] args) throws IOException
		{
		new Gtf2Xml().instanceMainWithExit(args);
		}
}
