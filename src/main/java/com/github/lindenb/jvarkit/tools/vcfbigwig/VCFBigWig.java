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

*/
package com.github.lindenb.jvarkit.tools.vcfbigwig;

import java.io.Closeable;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;
import javax.xml.bind.annotation.XmlType;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.events.XMLEvent;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.math.stats.Percentile;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFIterator;
/*
BEGIN_DOC

## XML definition

multiple BigWig files can be specified using a XML file.

* Root is `<registry>`
* under `<registry>` is a set of `<bigwig>' elements.
* under `<bigwig>` contains the `<uri>'(required) , `<tag>` and `<description>`


## Example

```bash
 java -jar dist/vcfbigwig.jar \
 	-T GERP \
 	-B gerp.bw input.vcf.gz 
	
##INFO=<ID=GERP,Number=1,Type=Float,Description="Values from bigwig file: com.github.lindenb.jvarkit.tools.vcfbigwig.VCFBigWig BIGWIG=gerp.bw TAG=GERP IN=input.vcf.gz    VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO(...)
A	33926	.	G	A	182	.	GERP=-6.35(...)
A	45365	.	A	G	222	.	GERP=-3.55(...)
```


END_DOC
*/
@Program(name="vcfbigwig",
	description="Annotate a VCF with values from a bigwig file",
	keywords={"vcf","wig","wiggle","bigwig"}
	)
public class VCFBigWig extends Launcher
	{
	private static final Logger LOG = Logger.build(VCFBigWig.class).make();

	/** wraps a BigWigIterator */
	private static class WigItemIterator
		extends AbstractIterator<WigItem>
		{
		private final BigWigIterator delegate;
		WigItemIterator(final BigWigIterator delegate ){
			this.delegate  = delegate;
			}
		@Override
		protected WigItem advance() {
			if(this.delegate==null) return null;
			if(!this.delegate.hasNext()) return null;
			return this.delegate.next();
			}
		}
	
	/** describe a BigWig Resource */
	private static class BigWigResource
		implements Closeable
		{
		private String tag;
		private String biwWigFile;
		private String description;
		private BBFileReader bbFileReader=null;
		private ContigNameConverter contigNameConverter = null;
		private final Set<String> userContigsNotFound = new HashSet<>();
		
	
		public String getToken() {
			if(StringUtil.isBlank(this.tag))
				{
				this.tag=this.biwWigFile;
				int i=this.tag.lastIndexOf(File.separator);
				if(i!=-1) this.tag=this.tag.substring(i+1);
				i=this.tag.indexOf('.');
				if(i!=-1) this.tag=this.tag.substring(0,i);
				if(StringUtil.isBlank(this.tag)) throw new JvarkitException.UserError("Bad TAG for "+this.biwWigFile);
				LOG.info("setting tag to "+this.tag+" for "+this.biwWigFile);
				}
			return this.tag;
			}
		
		public String getPath() { return this.biwWigFile;}
		
		public String getDescription() { return
				StringUtil.isBlank(this.description)?
					getPath():
					this.description
					;}
		
		
		public BigWigResource open()
			{
			try {
				this.bbFileReader= new BBFileReader(this.biwWigFile);
				}
			catch(final IOException err)
				{
				throw new RuntimeIOException("Cannot open "+this.biwWigFile,err);
				}
			if(!this.bbFileReader.isBigWigFile())
				{
				this.bbFileReader=null;
				throw new RuntimeIOException(this.biwWigFile+" is not a bigWIG file. ("+this.getToken()+")");
				}
			this.contigNameConverter = ContigNameConverter.fromContigSet(new HashSet<>(this.bbFileReader.getChromosomeNames()));
			return this;
			}
		
		public Iterator<WigItem> iterator(final Locatable locatable,boolean contained)
			{
			return new WigItemIterator(this.bbFileReader.getBigWigIterator(
					locatable.getContig(),
					locatable.getStart()-1,
					locatable.getContig(),
					locatable.getEnd(),
					contained
					));
			}
		
		@Override
		public void close() {
			try
				{
				if(this.bbFileReader!=null)
					{
					CloserUtil.close(this.bbFileReader.getBBFis());
					}
				CloserUtil.close(this.bbFileReader);
				this.bbFileReader=null;
				}
			catch(final Exception err)
				{
				LOG.error(err);
				}
			}
		}

	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@ParametersDelegate
	private CtxWriterFactory component = new CtxWriterFactory();

	@XmlType(name="vcfbigwig")
	@XmlRootElement(name="vcfbigwig")
	@XmlAccessorType(XmlAccessType.FIELD)
	public static class CtxWriterFactory 
	implements VariantContextWriterFactory
		{
		@XmlElement(name="bigwig")
		@Parameter(names={"-B","--bigwig"},description=
				"Path to the bigwig file. "
				+ "[20180122] If the path ends with '.xml' it is interpretted as a XML file containing describing a set of BigWig resources; See online doc.",required=true)
		private String userBigWigFileUri = null;
	
		@XmlElement(name="tag")
		@Parameter(names={"-T","--tag","-tag"},description="Name of the INFO tag. default: name of the bigwig")
		private String userVcfTag = null;
	
		@Parameter(names={"-C","--contained"},description="Specifies wig values must be contained by region. if false: return any intersecting region values")
		private boolean contained = false;
	
		@XmlElement(name="aggregate")
		@Parameter(names={"-a","--aggregate"},description="How to aggregate overlapping values: 'avg' average; 'median': median, 'first': use first, 'all' : print all the data")
		private AggregateMethod aggregateMethod  = AggregateMethod.avg;
	
		@XmlTransient
		@Parameter(names={"-t","--transform"},description="Deprecated",hidden=true)
		private String _convertChrName = null;
	
		
		private final List<BigWigResource> bigwigResources = new ArrayList<>();

		
		private class CtxWriter extends DelegateVariantContextWriter
			{
			private final AggregateMethod aggregateMethod;
			private final List<Float> values=new ArrayList<Float>();

			
			CtxWriter(final VariantContextWriter delegate) {
				super(delegate);
				this.aggregateMethod = CtxWriterFactory.this.aggregateMethod;
				}
			
			
			@Override
			public void writeHeader(final VCFHeader header) {
					
				final VCFHeader h2=new VCFHeader(header);
				
				for(final BigWigResource rsrc: CtxWriterFactory.this.bigwigResources) {
					
					if(h2.getInfoHeaderLine(rsrc.getToken())!=null)
						{
						throw new JvarkitException.DuplicateVcfHeaderInfo(h2,rsrc.getToken());
						}
					
					if(this.aggregateMethod.equals(AggregateMethod.all))
						{
						h2.addMetaDataLine(new VCFInfoHeaderLine(
								rsrc.getToken(),
								VCFHeaderLineCount.UNBOUNDED,
								VCFHeaderLineType.Float,
								"Values from bigwig file: "+rsrc.getPath()+". "+rsrc.getDescription()
								));
						}
					else
						{
						h2.addMetaDataLine(new VCFInfoHeaderLine(
								rsrc.getToken(),1,
								VCFHeaderLineType.Float,
								"Values from bigwig file: "+rsrc.getPath()+". "+rsrc.getDescription()
								));
						}
					}
				
				super.writeHeader(h2);
				}

			@Override
			public void add(final VariantContext ctx) {
				VariantContextBuilder vcb = null;
				
				
				for(final BigWigResource rsrc:CtxWriterFactory.this.bigwigResources) {
					this.values.clear();
					final String variantChrom=  rsrc.contigNameConverter.apply(ctx.getContig());
					
					if(StringUtil.isBlank(variantChrom)) {
						if(!rsrc.userContigsNotFound.contains(ctx.getContig()))
							{
							rsrc.userContigsNotFound.add(ctx.getContig());
							LOG.warn("Bigwig file \""+rsrc.getPath()+"\" doesn't contains contig "+ variantChrom+"/"+ctx.getContig());
							}
						continue;
						}
					
					
					final Iterator<WigItem> iter=rsrc.iterator(
							new Interval(variantChrom,ctx.getStart(),ctx.getEnd()),
							CtxWriterFactory.this.contained
							);
					while(iter!=null && iter.hasNext())
						{
						final WigItem item=iter.next();
						final float v=item.getWigValue();
						this.values.add(v);
						if(this.aggregateMethod.equals(AggregateMethod.first)) break;
						}
					
					if(this.values.isEmpty())
						{
						continue;
						}
					if(vcb==null) vcb=new VariantContextBuilder(ctx);
	
					switch(this.aggregateMethod)
						{
						case all:
							vcb.attribute(rsrc.getToken(),this.values);
							break;
						case avg:
							vcb.attribute(rsrc.getToken(),
									(float)Percentile.average().evaluate(values.stream().mapToDouble(V->V.doubleValue()).toArray()));
							break;
						case first:
							vcb.attribute(rsrc.getToken(),values.get(0));
							break;
						case median:
							vcb.attribute(rsrc.getToken(),
									(float)Percentile.median().evaluate(values.stream().mapToDouble(V->V.doubleValue()).toArray()));
							break;
						default: throw new IllegalStateException();
						}
					}
				if(vcb==null)
					{
					super.add(ctx);
					}
				else
					{
					super.add(vcb.make());
					}
				
				
				}
			
			@Override
			public void close() {
				for(final BigWigResource rsrc:CtxWriterFactory.this.bigwigResources)
					{
					if(!rsrc.userContigsNotFound.isEmpty())
						{
						LOG.warn("\""+rsrc.getPath()+
								"\": Contigs not found :"+
								String.join(" ", rsrc.userContigsNotFound));
						}
					rsrc.close();
					}
				super.close();
				}
			
			}
		
		@Override
		public int initialize() {
			FileReader fr =null;
			try
				{
				if(StringUtil.isBlank(this.userBigWigFileUri))
					{
					LOG.info("Undefined BigWig file ");
					return -1;
					}
				
				if(this.userBigWigFileUri.endsWith(".xml"))
					{
					XMLInputFactory xif=XMLInputFactory.newFactory();
					fr=new FileReader(new File(this.userBigWigFileUri));
					XMLEventReader r = xif.createXMLEventReader(fr);
					XMLEvent evt=r.nextTag();
					if(evt==null || !evt.isStartElement() ||
							!evt.asStartElement().getName().getLocalPart().equals("registry"))
						{
						LOG.error("Root of "+this.userBigWigFileUri+" is not <registry> "+evt);
						}
					while(r.hasNext())
						{
						evt=r.nextEvent();
						if(evt.isEndElement()) break;
						if(!evt.isStartElement()) continue;
						if(!evt.asStartElement().getName().getLocalPart().equals("bigwig")) continue;
						final BigWigResource rsrc = new BigWigResource();
						while(r.hasNext())
							{
							final XMLEvent evt2=r.nextEvent();
							if(evt2.isEndElement()) break;
							if(!evt2.isStartElement()) continue;
							final String localName= evt2.asStartElement().getName().getLocalPart();
							if(localName.equals("uri") && rsrc.biwWigFile==null)
								{
								rsrc.biwWigFile=r.getElementText().trim();
								}
							else if(localName.equals("tag") && rsrc.tag==null)
								{
								rsrc.tag =r.getElementText().trim();
								}
							else if(localName.equals("description") && rsrc.description==null)
								{
								rsrc.description=r.getElementText().trim();
								}
							else
								{
								LOG.error("BAD XML element "+evt2.getLocation()+" "+localName);
								rsrc.close();
								return -1;
								}	
							}
						if(StringUtil.isBlank(rsrc.biwWigFile)) {
								LOG.error("undefined uri in "+evt.getLocation());
								rsrc.close();
								return -1;
							}
						if(this.bigwigResources.stream().anyMatch(R->R.getToken().equals(rsrc.getToken())))
							{
							LOG.error("BigWig Named "+rsrc.getToken()+" defined twice");
							rsrc.close();
							return -1;
							}
						this.bigwigResources.add(rsrc);
						}
					
					r.close();
					fr.close();fr=null;
					}
				else
					{
					final BigWigResource rsrc = new BigWigResource();
					rsrc.biwWigFile = this.userBigWigFileUri;
					rsrc.tag = this.userVcfTag;
					rsrc.description=this.userBigWigFileUri;
					this.bigwigResources.add(rsrc);
					}
				
				this.bigwigResources.stream().forEach(BB->{
					BB.open();
					});
				return 0;
				}
			catch(final Exception err)
				{
				LOG.error(err);
				return -1;
				}
			finally
				{
				CloserUtil.close(fr);
				}
			}
		
		@Override
		public VariantContextWriter open(VariantContextWriter delegate) {
			return new CtxWriter(delegate);
			}

		@Override
		public void close() throws IOException {
			CloserUtil.close(this.bigwigResources);
			VariantContextWriterFactory.super.close();
			}
		}
	
	
	private enum AggregateMethod
		{
		avg,median,first,all
		}
	public VCFBigWig()
		{
		}
	
	
	@Override
	protected int doVcfToVcf(final String inputName, final VCFIterator r, final VariantContextWriter delegate) {
		final VariantContextWriter w = this.component.open(delegate);
		w.writeHeader(r.getHeader());
		final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(r.getHeader()).logger(LOG);
		while(r.hasNext())
			{
			w.add(progress.watch(r.next()));
			
			// JVM crash sometimes ? suspect there is a memory leak ?
			if(progress.getCount()%1000L==0)
				{
				System.gc();
				}
			}
		progress.finish();
		w.close();
		return 0;
		}
	
	@Override
	public int doWork(final List<String> args) {
		try 
			{
			if(this.component.initialize()!=0) {
				return -1;
				}
			return doVcfToVcf(args,this.outputFile);
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.component);
			}
		}
	
	public static void main(final String[] args) throws IOException
		{
		new VCFBigWig().instanceMainWithExit(args);
		}
}
