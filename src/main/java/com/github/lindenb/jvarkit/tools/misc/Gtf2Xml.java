/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.readers.LineIterator;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFCodec;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFLine;
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
	description="Convert GTF/GFF to XML",
	keywords={"xml","gtf","gff","gff3"})
public class Gtf2Xml extends Launcher{
	private static final Logger LOG = Logger.build(FixVCF.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-f","--features"},description="Don't record features types.")
	private boolean disable_feature_type=false;
	@Parameter(names={"-s","--sources"},description="Don't record sources")
	private boolean disable_sources=false;
	@Parameter(names={"-d","--dict"},description="Don't record contigs")
	private boolean disable_dict=false;
	@Parameter(names={"-a","--attributes"},description="Don't record attribute types.")
	private boolean disable_att_keys=false;
	@ParametersDelegate
	private GTFCodec.FormatChooser formatChooser = new  GTFCodec.FormatChooser();
	

	
	
	private final Map<String,Long> seqdict=new LinkedHashMap<>();
	private final Set<String> att_keys=new HashSet<>();
	private final Set<String> sources=new HashSet<>();
	private final Set<String> types=new HashSet<>();
	
	private void write(final XMLStreamWriter w,final GTFLine line) throws XMLStreamException,IOException
		{
		w.writeStartElement("feature");
		
		w.writeAttribute("chrom", line.getContig());
		w.writeAttribute("start",String.valueOf(line.getStart()));
		w.writeAttribute("end",String.valueOf(line.getEnd()));
		if(!this.disable_dict) {
			Long contifLength  = this.seqdict.get(line.getContig());
			if(contifLength==null) contifLength=0L;
			this.seqdict.put(line.getContig(),Math.max(line.getEnd(), contifLength));
			}
		if(line.getScore()!=null)  {
			w.writeAttribute("score",String.valueOf(line.getScore()));
		}
		if(line.getStrand()!=GTFLine.NO_STRAND) {
			w.writeAttribute("strand",String.valueOf(line.getStrand()));
		}
		if(!StringUtil.isBlank(line.getSource()))  {
			if(!this.disable_sources) this.sources.add(line.getSource());
			w.writeAttribute("source",line.getSource());
		}
		if(!StringUtil.isBlank(line.getType()))  {
			if(!disable_feature_type) this.types.add(line.getType());
			w.writeAttribute("type",line.getType());
		}
		if(line.getPhase()!=GTFLine.NO_PHASE) {
			w.writeAttribute("phase",String.valueOf(line.getPhase()));
			}
		
		w.writeStartElement("attributes");
		for(final Iterator<Map.Entry<String, String>> kvr=line.getAttributeIterator();
				kvr.hasNext();
				)
			{
			final Map.Entry<String, String> kv = kvr.next();
			w.writeStartElement("attribute");
			w.writeAttribute("key", kv.getKey());
			if(!this.disable_att_keys) this.att_keys.add(kv.getKey());
			w.writeCharacters(kv.getValue());
			w.writeEndElement();
			}
		w.writeEndElement();
			
		w.writeEndElement();
		w.writeCharacters("\n");
		}
	
	@Override
	public int doWork(final List<String> args) {

		LineIterator r=null;
		XMLStreamWriter w=null;
		FileWriter fw=null;
		try {
			String inputName=oneFileOrNull(args);
			r = (StringUtil.isBlank(inputName)?
					IOUtils.openStreamForLineIterator(stdin()):
					IOUtils.openURIForLineIterator(inputName)
					);
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			if(this.outputFile==null)
				{
				w = xof.createXMLStreamWriter(stdout(), "UTF-8");
				}
			else
				{
				w = xof.createXMLStreamWriter((fw=new FileWriter(this.outputFile)));
				}
			final GTFCodec codec = this.formatChooser.makeCodec();
			w.writeStartDocument("UTF-8","1.0");
			w.writeStartElement("gtf");
			
			final GTFCodec.GTFHeader header = codec.readActualHeader(r);
			for(final String headerLine:header.getLines()) {
				if(!headerLine.startsWith("#!")) continue;
				final int ws = headerLine.indexOf(' ');
				if(ws==-1) continue; //??
				w.writeAttribute(
						headerLine.substring(2, ws),
						headerLine.substring(ws+1).trim());
				
				}
		
			while(r.hasNext())
				{
				final String line=r.next();
				GTFLine gtfline = codec.decode(line);
				if(gtfline==null) continue;
				write(w,gtfline);
				}
			
			if(!this.disable_att_keys) {
				w.writeStartElement("attributes");
				for(String k : this.att_keys)
					{
					w.writeStartElement("attribute");
					w.writeCharacters(k);
					w.writeEndElement();
					w.writeCharacters("\n");
					}
				w.writeEndElement();
				w.writeCharacters("\n");
				}
			
			
			if(!this.disable_feature_type) {
				w.writeStartElement("types");
				for(String k : this.types)
					{
					w.writeStartElement("type");
					w.writeCharacters(k);
					w.writeEndElement();
					w.writeCharacters("\n");
					}
				w.writeEndElement();
				w.writeCharacters("\n");
				}
			
			if(!this.disable_sources) {
				w.writeStartElement("sources");
				for(final String k : this.sources)
					{
					w.writeStartElement("source");
					w.writeCharacters(k);
					w.writeEndElement();
					w.writeCharacters("\n");
					}
				w.writeEndElement();
				w.writeCharacters("\n");
				}

			
			if(!this.disable_dict) {
				w.writeStartElement("dict");
				for(final String k : this.seqdict.keySet())
					{
					w.writeEmptyElement("chrom");
					w.writeAttribute("name",k);
					w.writeAttribute("length", String.valueOf(this.seqdict.get(k)));
					w.writeCharacters("\n");
					}
				w.writeEndElement();
				w.writeCharacters("\n");
				}

			
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
