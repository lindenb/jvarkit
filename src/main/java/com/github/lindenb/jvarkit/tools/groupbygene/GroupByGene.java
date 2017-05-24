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
package com.github.lindenb.jvarkit.tools.groupbygene;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;

/**
BEGIN_DOC

## Motivation

Group VCF data by gene/transcript. By default it tries to use data from VEP and SnpEff

## Example

### Delimited output

```
$ curl -s -k "https://raw.github.com/arq5x/gemini/master/test/test4.vep.snpeff.vcf" |\
java -jar dist/groupbygene.jar |\
head | column  -t

#chrom  min.POS    max.POS    gene.name  gene.type         samples.affected  count.variations  M10475  M10478  M10500  M128215
chr10   52004315   52004315   ASAH2      snpeff-gene-name  2                 1                 0       0       1       1
chr10   52004315   52004315   ASAH2      vep-gene-name     2                 1                 0       0       1       1
chr10   52497529   52497529   ASAH2B     snpeff-gene-name  2                 1                 0       1       1       0
chr10   52497529   52497529   ASAH2B     vep-gene-name     2                 1                 0       1       1       0
chr10   48003992   48003992   ASAH2C     snpeff-gene-name  3                 1                 1       1       1       0
chr10   48003992   48003992   ASAH2C     vep-gene-name     3                 1                 1       1       1       0
chr10   126678092  126678092  CTBP2      snpeff-gene-name  1                 1                 0       0       0       1
chr10   126678092  126678092  CTBP2      vep-gene-name     1                 1                 0       0       0       1
chr10   135336656  135369532  CYP2E1     snpeff-gene-name  3                 2                 0       2       1       1
```

## XML output

```
$ curl -s -k "https://raw.github.com/arq5x/gemini/master/test/test4.vep.snpeff.vcf" |\
java -jar dist/groupbygene.jar -X |\
xmllint --> --format -<!-- 
```


```
<!-- <?xml version="1.0" encoding="UTF-8"?>
<genes>
  <samples count="4">
    <sample>M10475</sample>
    <sample>M10478</sample>
    <sample>M10500</sample>
    <sample>M128215</sample>
  </samples>
  <gene name="ASAH2" type="snpeff-gene-name" chrom="chr10" min.POS="52004315" max.POS="52004315" affected="2" variations="1">
    <sample name="M10500" count="1">
      <genotype pos="52004315" ref="T" A1="C" A2="C"/>
    </sample>
    <sample name="M128215" count="1">
      <genotype pos="52004315" ref="T" A1="C" A2="C"/>
    </sample>
  </gene>
  <gene name="ASAH2" type="vep-gene-name" chrom="chr10" min.POS="52004315" max.POS="52004315" affected="2" variations="1">
    <sample name="M10500" count="1">
(...)
    <sample name="M10475" count="1">
      <genotype pos="72057435" ref="C" A1="C" A2="T"/>
    </sample>
  </gene>
  <gene name="ENST00000572003" type="vep-ensembl-transcript-name" chrom="chr16" min.POS="72057435" max.POS="72057435" affected="1" variations="1">
    <sample name="M10475" count="1">
      <genotype pos="72057435" ref="C" A1="C" A2="T"/>
    </sample>
  </gene>
  <gene name="ENST00000572887" type="vep-ensembl-transcript-name" chrom="chr16" min.POS="72057435" max.POS="72057435" affected="1" variations="1">
    <sample name="M10475" count="1">
      <genotype pos="72057435" ref="C" A1="C" A2="T"/>
    </sample>
  </gene>
  <gene name="ENST00000573843" type="vep-ensembl-transcript-name" chrom="chr16" min.POS="72057435" max.POS="72057435" affected="1" variations="1">
    <sample name="M10475" count="1">
      <genotype pos="72057435" ref="C" A1="C" A2="T"/>
    </sample>
  </gene>
  <gene name="ENST00000573922" type="vep-ensembl-transcript-name" chrom="chr16" min.POS="72057435" max.POS="72057435" affected="1" variations="1">
    <sample name="M10475" count="1">
      <genotype pos="72057435" ref="C" A1="C" A2="T"/>
    </sample>
  </gene>
  <gene name="ENST00000574309" type="vep-ensembl-transcript-name" chrom="chr16" min.POS="72057435" max.POS="72057435" affected="1" variations="1">
    <sample name="M10475" count="1">
      <genotype pos="72057435" ref="C" A1="C" A2="T"/>
    </sample>
  </gene>
</genes>
```

END_DOC
 */
@Program(name="groupbygene",keywords={"vcf","gene"},description="Group VCF data by gene/transcript. By default it uses data from VEP , SnpEff")
public class GroupByGene
	extends Launcher
	{
	private static final Logger LOG = Logger.build(GroupByGene.class).make();

	@Parameter(names={"-X","--xml"},description="XML output")
	private boolean xml_output = false;
	@Parameter(names={"--filtered"},description="ignore FILTERED variants")
	private boolean ignore_filtered = false;
	@Parameter(names={"-T","--tag"},description="add Tag in INFO field containing the name of the genes.")
	private Set<String> user_gene_tags = new HashSet<>();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outFile=null;
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	
	private Set<String> sampleNames = new TreeSet<String>();
	private SortingCollection<Call> sortingCollection = null;
	
	private static class GeneName
		{
		String name;
		String type;
		GeneName(String name,String type)
			{
			this.name=name;
			this.type=type;
			}
		@Override
		public int hashCode()
			{
			final int prime = 31;
			int result = 1;
			result = prime * result +  name.hashCode();
			result = prime * result +  type.hashCode();
			return result;
			}
		@Override
		public boolean equals(Object o)
			{
			if (this == o) return true;
			if (o == null) return false;
			if (getClass() != o.getClass()) return false;
			GeneName g=(GeneName)o;
			return name.equals(g.name) && type.equals(g.type);
			}
		@Override
		public String toString() {
			return  name+"("+type+")";
			}
		
		}
	
	private static class Call
		{
		String chrom;
		int pos;
		String ref;
		GeneName gene;
		String sample;
		String posId;
		String a1=null;
		String a2=null;
		}
	
	private static class CallCodec
		extends AbstractDataCodec<Call>
		{
		@Override
		public void encode(DataOutputStream dos, Call c)
				throws IOException
			{
			dos.writeUTF(c.chrom);
			dos.writeInt(c.pos);
			dos.writeUTF(c.ref);
			dos.writeUTF(c.gene.name);
			dos.writeUTF(c.gene.type);
			
			dos.writeUTF(c.sample);
			dos.writeUTF(c.posId);
			dos.writeUTF(c.a1);
			dos.writeUTF(c.a2);
			}
		
		@Override
		public Call decode(DataInputStream dis) throws IOException
			{
			final Call c= new Call();
			try {
				c.chrom=dis.readUTF();
			} catch (Exception e) {
				return null;
				}
			c.pos=dis.readInt();
			c.ref=dis.readUTF();
			String gName=dis.readUTF();
			String gType=dis.readUTF();
			c.gene=new GeneName(gName, gType);
			c.sample=dis.readUTF();
			c.posId=dis.readUTF();
			c.a1=dis.readUTF();
			c.a2=dis.readUTF();
			return c;
			}
		@Override
		public CallCodec clone() {
			return new CallCodec();
			}
		}
	
	private static class CallCmp
		implements Comparator<Call>
		{
		@Override
		public int compare(final Call o1, final Call o2)
			{
			int i=o1.chrom.compareTo(o2.chrom);
			if(i!=0) return i;
			 i= o1.gene.name.compareTo(o2.gene.name);
			 if(i!=0) return i;
			 i= o1.gene.type.compareTo(o2.gene.type);
			 return i;
			}
		}
	
	private SnpEffPredictionParser snpEffPredictionParser=null;
	private VepPredictionParser vepPredictionParser=null;
	private AnnPredictionParser annPredictionParser=null;
	
	public GroupByGene()
		{
		}
	
	private Set<GeneName> getGenes(final VariantContext ctx)
		{
		final Set<GeneName> set=new HashSet<GeneName>();
		for(VepPredictionParser.VepPrediction pred: this.vepPredictionParser.getPredictions(ctx))
			{
			String s=pred.getGeneName();
			if(s!=null)  set.add(new GeneName(s,"vep-gene-name"));
			s=pred.getEnsemblGene();
			if(s!=null)  set.add(new GeneName(s,"vep-ensembl-gene-name"));
			s=pred.getSymbol();
			if(s!=null)  set.add(new GeneName(s,"vep-symbol"));
			s=pred.getHGNC();
			if(s!=null)  set.add(new GeneName(s,"vep-hgnc"));
			s=pred.getHgncId();
			if(s!=null)  set.add(new GeneName(s,"vep-hgnc-id"));
			s=pred.getGene();
			if(s!=null)  set.add(new GeneName(s,"vep-gene"));
			s=pred.getRefSeq();
			if(s!=null)  set.add(new GeneName(s,"vep-refseq"));
			}
		for(SnpEffPredictionParser.SnpEffPrediction pred: this.snpEffPredictionParser.getPredictions(ctx))
			{
			String s=pred.getGeneName();
			if(s!=null)  set.add(new GeneName(s,"snpeff-gene-name"));
			s=pred.getEnsemblTranscript();
			if(s!=null)  set.add(new GeneName(s,"snpeff-ensembl-transcript-name"));
			}
		for(final AnnPredictionParser.AnnPrediction pred:this.annPredictionParser.getPredictions(ctx)) {
			String s=pred.getGeneId();
			if(s!=null)  set.add(new GeneName(s,"ann-gene-id"));
			s=pred.getGeneName();
			if(s!=null)  set.add(new GeneName(s,"ann-gene-name"));
			s=pred.getFeatureId();
			if(s!=null)  set.add(new GeneName(s,"ann-feature-id"));
			}
		
		for(String user_gene_tag:user_gene_tags)
			{
			if(user_gene_tag.isEmpty()) continue;
			if(user_gene_tag.equals(".")) continue;
			Object o=ctx.getAttribute(user_gene_tag);
			Set<String> tag=new HashSet<String>();
			if(o==null) 
				{
				//
				}
			else if(o.getClass().isArray())
				{
				for(Object o2:(Object[])o) tag.add(String.valueOf(o2));
				}
			else if(o instanceof java.util.Collection)
				{
				for(Object o2:(java.util.Collection<?>)o) tag.add(String.valueOf(o2));
				}
			else
				{
				tag.add(o.toString());
				}
			for(String t:tag)
				{
				if(t.isEmpty()) continue;
				set.add(new GeneName(t,"user:"+user_gene_tag));
				}
			}
		Iterator<GeneName> iter=set.iterator();
		while(iter.hasNext())
			{
			if(iter.next().name.isEmpty()) iter.remove();
			}
		return set;
		}
	
	public void addUserGeneTag(String t) {
		this.user_gene_tags.add(t);
		}
	
	/** public for knime */
	public void dump() throws IOException,XMLStreamException
		{
		this.sortingCollection.doneAdding();
		
		PrintStream pw = openFileOrStdoutAsPrintStream(this.outFile);
		
		XMLStreamWriter w=null;
		if(xml_output)
			{
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			w= xof.createXMLStreamWriter(pw, "UTF-8");
			w.writeStartDocument("UTF-8","1.0");
			w.writeStartElement("genes");
			w.writeComment("Cmd line: "+this.getProgramCommandLine());
			w.writeComment("Version "+getVersion());
			w.writeStartElement("samples");
			w.writeAttribute("count", String.valueOf(sampleNames.size()));
			for(String sample:this.sampleNames)
				{
				w.writeStartElement("sample");
				w.writeCharacters(sample);
				w.writeEndElement();
				}
			w.writeEndElement();
			w.writeCharacters("\n");
			}
		else
			{
			pw.print("#chrom");
			pw.print('\t');
			pw.print("min.POS");
			pw.print('\t');
			pw.print("max.POS");
			pw.print('\t');
			pw.print("gene.name");
			pw.print('\t');
			pw.print("gene.type");
			pw.print('\t');
			pw.print("samples.affected");
			pw.print('\t');
			pw.print("count.variations");
			for(String sample:this.sampleNames)
				{
				pw.print('\t');
				pw.print(sample);
				}
			pw.println();
			}
		
		final CallCmp cmp=new CallCmp();
		List<Call> row=new ArrayList<Call>();
		CloseableIterator<Call> iter=sortingCollection.iterator();
		for(;;)
			{
			Call curr=null;
			if(iter.hasNext()) curr=iter.next();
			if(curr==null || (!row.isEmpty() && cmp.compare(curr, row.get(0))!=0))
				{
				if(!row.isEmpty())
					{
					int minPos=Integer.MAX_VALUE;
					int maxPos=Integer.MIN_VALUE;
					Set<String> affected=new HashSet<String>();
					Set<String> distinctMut=new HashSet<String>();
					Counter<String> sample2count=new Counter<String>();
					for(Call c:row)
						{
						minPos=Math.min(minPos, c.pos);
						maxPos=Math.max(maxPos, c.pos);
						sample2count.incr(c.sample);
						affected.add(c.sample);
						distinctMut.add(c.chrom+":"+c.pos+":"+c.ref);
						}
					Call first=row.get(0);
					if(w!=null)
						{
						w.writeStartElement("gene");
						w.writeAttribute("name", first.gene.name);
						w.writeAttribute("type", first.gene.type);
						
						w.writeAttribute("chrom", first.chrom);
						w.writeAttribute("min.POS",String.valueOf(minPos));
						w.writeAttribute("max.POS",String.valueOf(maxPos));
						w.writeAttribute("affected",String.valueOf(affected.size()));
						w.writeAttribute("variations",String.valueOf(distinctMut.size()));
						
						for(String sample:this.sampleNames)
							{
							if(sample2count.count(sample)==0L) continue;
							w.writeStartElement("sample");
							w.writeAttribute("name",sample);
							
							w.writeAttribute("count",String.valueOf(sample2count.count(sample)));
							for(Call c:row)
								{
								if(!c.sample.equals(sample)) continue;
								w.writeEmptyElement("genotype");
								w.writeAttribute("pos", String.valueOf(c.pos));
								w.writeAttribute("ref", c.ref);
								w.writeAttribute("A1", c.a1);
								w.writeAttribute("A2", c.a2);
								}
							w.writeEndElement();							
							}
						
						
						w.writeEndElement();//gene
						w.writeCharacters("\n");
						}
					else
						{
						pw.print(first.chrom);
						pw.print('\t');
						pw.print(minPos);
						pw.print('\t');
						pw.print(maxPos);
						pw.print('\t');
						pw.print(first.gene.name);
						pw.print('\t');
						pw.print(first.gene.type);
						pw.print('\t');
						pw.print(affected.size());
						pw.print('\t');
						pw.print(distinctMut.size());
						for(String sample:this.sampleNames)
							{
							pw.print('\t');
							pw.print(sample2count.count(sample));
							}
						pw.println();
						if(pw.checkError()) break;
						}
					}
				if(curr==null) break;
				row.clear();
				}
			row.add(curr);
			}
		iter.close();
		if(w!=null)
			{
			w.writeEndElement();
			w.writeEndDocument();
			w.flush();
			w.close();
			}
		pw.flush();
		if(this.outFile!=null) pw.close();
		}
	
	private void read(final InputStream in) throws IOException
		{
		final VcfIterator iter = VCFUtils.createVcfIteratorFromInputStream(in);
		this.readVcf(iter);
		CloserUtil.close(iter);
		}
	
	/** public for knime */
	public void readVcf(final VcfIterator iter) 
		{
		final VCFHeader header=(VCFHeader)iter.getHeader();
		if(header.getSampleNamesInOrder()!=null)
			{
			this.sampleNames.addAll(header.getSampleNamesInOrder());
			}
		this.snpEffPredictionParser=new SnpEffPredictionParserFactory(header).get();
		this.vepPredictionParser=new VepPredictionParserFactory(header).get();
		this.annPredictionParser = new AnnPredictionParserFactory(header).get();
		SAMSequenceDictionary dict=header.getSequenceDictionary();
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
		while(iter.hasNext())
			{
			final VariantContext ctx= progress.watch(iter.next());
			if(ignore_filtered && ctx.isFiltered()) continue;
			
			for(final GeneName g:getGenes(ctx))
				{
				for(final Genotype genotype:ctx.getGenotypes())
					{
					if(!genotype.isAvailable()) continue;
					if(!genotype.isCalled()) continue;
					if(genotype.isNoCall()) continue;
					if(genotype.isHomRef()) continue;
					final List<Allele> L=genotype.getAlleles();
					if(L==null || L.isEmpty()) continue;
					

					final Call c=new Call();
					c.chrom=ctx.getContig();
					c.pos=ctx.getStart();
					c.ref=ctx.getReference().getDisplayString();
					c.gene=g;
					c.sample=genotype.getSampleName();
					c.posId=(ctx.getContig()+":"+ctx.getStart()+":"+ctx.getReference().getDisplayString()).toLowerCase();
					if(L.size()==1)
						{
						c.a1=genotype.getAllele(0).getDisplayString().toUpperCase();
						c.a2=c.a1;
						}
					else if(L.size()==2)
						{
						c.a1=genotype.getAllele(0).getDisplayString().toUpperCase();
						c.a2=genotype.getAllele(1).getDisplayString().toUpperCase();
						if(c.a1.compareTo(c.a2)>0)
							{
							String tmp=c.a1;
							c.a1=c.a2;
							c.a2=tmp;
							}
						}
					else
						{
						CloserUtil.close(iter);
						throw new RuntimeException("cannot handle multi-ploidy "+ctx);
						}

					this.sortingCollection.add(c);
					}
				}	
			}
		CloserUtil.close(iter);
		}
	
	@Override
	public int doWork(List<String> args) {
		try
			{
			this.sortingCollection=SortingCollection.newInstance(
					Call.class,
					new CallCodec(),
					new CallCmp(),
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpDirectories()
					);
			this.sortingCollection.setDestructiveIteration(true);
			this.sampleNames.clear();
			
			if(args.isEmpty())
				{
				LOG.info("Reading from stdin");
				read(stdin());
				}
			else
				{
				for(final String filename:args)
					{
					LOG.info("Reading from "+filename);
					InputStream in=IOUtils.openURIForReading(filename);
					read(in);
					in.close();
					}
				}
	
			LOG.info("Done reading. Now printing results.");
	
			dump();
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			if(sortingCollection!=null) sortingCollection.cleanup();
			}
		}

	
	
	
	public static void main(String[] args)
		{
		new GroupByGene().instanceMainWithExit(args);
		}
	}
