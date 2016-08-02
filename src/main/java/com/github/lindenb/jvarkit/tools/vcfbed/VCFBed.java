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
package com.github.lindenb.jvarkit.tools.vcfbed;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import com.github.lindenb.jvarkit.util.bio.bed.IndexedBedReader;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


/**
 * 
 * VCFBed
 *
 */
public class VCFBed extends AbstractVCFBed
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VCFBed.class);
	private IntervalTreeMap<Set<BedLine>> intervalTreeMap=null;

	
	private IndexedBedReader bedReader =null;
	private Chunk parsedFormat=null;
	
	private static abstract class Chunk
		{
		public abstract String toString(BedLine tokens);
		public Chunk next=null;
		}
	
	private static class PlainChunk extends Chunk
		{
		final String s;
		PlainChunk(final String s){this.s=s;}
		public String toString(final BedLine tokens)
			{
			return s+(next==null?"":next.toString(tokens));
			}
		@Override
		public String toString() {
			return "plain-text:\""+this.s+"\""+(this.next==null?"":";"+this.next.toString());
			}
		}
	private static class ColChunk extends Chunk
		{
		final int index;
		ColChunk(final int index){ this.index=index;}
		public String toString(final BedLine tokens)
			{
			String s= tokens.get(index);
			if(s==null) s="";
			return s+(next==null?"":next.toString(tokens));
			}
		@Override
		public String toString() {
			return "column:${"+(index+1)+"}"+(this.next==null?"":";"+this.next.toString());
			}
		}

	
	private Chunk parseFormat(final String s)
		{
		if(s==null || s.isEmpty()) return null;
		if(s.startsWith("${"))
			{
			final int j=s.indexOf('}',2);
			if(j==-1) throw new IllegalArgumentException("bad format in \""+s+"\".");
			try
				{
				final int col=Integer.parseInt(s.substring(2, j).trim());
				if(col<1) throw new IllegalArgumentException();
				final ColChunk c=new ColChunk(col-1);
				c.next=parseFormat(s.substring(j+1));
				return c;
				}
			catch(final Exception err)
				{
				 throw new IllegalArgumentException("bad format in \""+s+"\".",err);
				}
			}
		else if(s.startsWith("$"))
			{
			int j=1;
			while(j<s.length() && Character.isDigit(s.charAt(j)))
				{
				++j;
				}
			int col=Integer.parseInt(s.substring(1, j).trim());
			if(col<1) throw new IllegalArgumentException();
			ColChunk c=new ColChunk(col-1);
			c.next=parseFormat(s.substring(j));
			return c;
			}
		int i=0;
		final StringBuilder sb=new StringBuilder();
		while(i< s.length() && s.charAt(i)!='$')
			{
			sb.append(s.charAt(i));
			i++;
			}
		final PlainChunk c=new PlainChunk(sb.toString());
		c.next=parseFormat(s.substring(i));
		return c;
		}
	
	
	@Override
	protected Collection<Throwable> doVcfToVcf(String inputName, VcfIterator r, VariantContextWriter w)
			throws IOException {
		final File srcbedfile = super.tabixFile==null?super.treeMapFile:super.tabixFile;
		final VCFHeader h2=new VCFHeader(r.getHeader());
		final VCFInfoHeaderLine infoHeader= 
				new VCFInfoHeaderLine(
						super.infoName,
						VCFHeaderLineCount.UNBOUNDED,
						VCFHeaderLineType.String,
						"metadata added from "+ srcbedfile+
						" . Format was "+super.formatPattern
						);
		
		final VCFFilterHeaderLine filterOverlap = 
				(this.filterOverlapStr==null || this.filterNoOverlapStr.trim().isEmpty()?null:
				new VCFFilterHeaderLine(this.filterOverlapStr, "Variant overlap with "+srcbedfile)	
				);
		
		final VCFFilterHeaderLine filterNoOverlap = 
				(this.filterNoOverlapStr==null || this.filterNoOverlapStr.trim().isEmpty()?null:
				new VCFFilterHeaderLine(this.filterNoOverlapStr, "Variant having NO overlap with "+srcbedfile)	
				);
		
		if(filterOverlap!=null) h2.addMetaDataLine(filterOverlap);
		if(filterNoOverlap!=null) h2.addMetaDataLine(filterNoOverlap);
		
		h2.addMetaDataLine(infoHeader);
		addMetaData(h2);
		
		final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(h2);
		w.writeHeader(h2);
		while(r.hasNext())
			{
			boolean found_overlap=false;
			final VariantContext ctx= progress.watch(r.next());
			final Set<String> annotations=new HashSet<String>();
			
			if(this.intervalTreeMap!=null) {
				for(final Set<BedLine> bedLines :this.intervalTreeMap.getOverlapping(new Interval(ctx.getContig(),ctx.getStart(),ctx.getEnd()))) {
					for(final BedLine bedLine:bedLines) {
					final String newannot=this.parsedFormat.toString(bedLine);
					found_overlap=true;
					if(!newannot.isEmpty())
						{
						annotations.add(VCFUtils.escapeInfoField(newannot));
						}
					}
				}
				
			}
			else
				{
				CloseableIterator<BedLine> iter = this.bedReader.iterator(
						ctx.getContig(),
						ctx.getStart()-1,
						ctx.getEnd()+1
						);
				while(iter.hasNext())
					{
					final BedLine bedLine = iter.next();
					
					if(!ctx.getContig().equals(bedLine.getContig())) continue;
					if(ctx.getStart() > bedLine.getEnd() ) continue;
					if(ctx.getEnd() < bedLine.getStart() ) continue;
	
					found_overlap=true;

					final String newannot=this.parsedFormat.toString(bedLine);
					if(!newannot.isEmpty())
						annotations.add(VCFUtils.escapeInfoField(newannot));
					}
				CloserUtil.close(iter);
				}
			
			final String filterToSet;
			if(found_overlap && filterOverlap!=null) {
				filterToSet = filterOverlap.getID();
				}
			else if(!found_overlap &&  filterNoOverlap!=null) {
				filterToSet = filterNoOverlap.getID();
				}
			else
				{
				filterToSet=null;
				}
			
			if(filterToSet==null && annotations.isEmpty())
				{
				w.add(ctx);
				continue;
				}
			final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
			if(!annotations.isEmpty()) {
				vcb.attribute(infoHeader.getID(), annotations.toArray());
				}
			if(filterToSet!=null) {
				vcb.filter(filterToSet);
				}
			w.add(vcb.make());
			if(w.checkError()) break;
			}
		progress.finish();
		return RETURN_OK;
		}
	
	/** reads a Bed file and convert it to a IntervalTreeMap<Bedline> */
	private htsjdk.samtools.util.IntervalTreeMap<Set<com.github.lindenb.jvarkit.util.bio.bed.BedLine>> readBedFileAsIntervalTreeMap(final java.io.File file) throws java.io.IOException
		{
		java.io.BufferedReader r=null;
		try
			{
			final  htsjdk.samtools.util.IntervalTreeMap<Set<com.github.lindenb.jvarkit.util.bio.bed.BedLine>> intervals = new
					 htsjdk.samtools.util.IntervalTreeMap<>();
			r=com.github.lindenb.jvarkit.io.IOUtils.openFileForBufferedReading(file);
			String line;
			final com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec codec = new com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec();
			while((line=r.readLine())!=null) 
				{
				if(line.startsWith("#") ||  com.github.lindenb.jvarkit.util.bio.bed.BedLine.isBedHeader(line) ||  line.isEmpty()) continue; 
				final com.github.lindenb.jvarkit.util.bio.bed.BedLine bl = codec.decode(line);
				if(bl==null || bl.getStart()>bl.getEnd()) continue;
				final htsjdk.samtools.util.Interval interval= bl.toInterval();
				Set<BedLine> set = intervals.get(interval);
				if(set==null )
					{
					set = new HashSet<>();
					intervals.put(interval,set); 
					}
				set.add(bl);
				}
			return intervals;
			}
		finally
			{
			htsjdk.samtools.util.CloserUtil.close(r);
			}
		}

	
	
	@Override
	public Collection<Throwable> initializeKnime() {
		try
			{
			if(super.tabixFile==null && super.treeMapFile==null)
				{
				return wrapException("Undefined tabix or memory file -"+OPTION_TABIXFILE+" -"+OPTION_TREEMAPFILE);
				}
			else if(super.tabixFile!=null && super.treeMapFile!=null)
				{
				return wrapException("You cannot use both options: -"+OPTION_TABIXFILE+" -"+OPTION_TREEMAPFILE);
				}
			else if( this.tabixFile!=null) {
				LOG.info("opening Bed "+this.tabixFile);
				this.bedReader= new IndexedBedReader(this.tabixFile);
				}
			else 
				{
				try {
					this.intervalTreeMap = this.readBedFileAsIntervalTreeMap(super.treeMapFile);
					LOG.info("Number of items in "+this.treeMapFile+" "+this.intervalTreeMap.size());
				}
				catch(final Exception err) {
					return wrapException(err);
				}
				}
			
			if(this.infoName==null || this.infoName.trim().isEmpty())
				{
				return wrapException("Undefined INFO name. -"+OPTION_INFONAME);
				}
			
			LOG.info("parsing "+this.formatPattern);
			this.parsedFormat=parseFormat(formatPattern);
			if(this.parsedFormat==null) this.parsedFormat=new PlainChunk("");
			LOG.info("format for "+this.formatPattern+" :"+this.parsedFormat);
			}
		catch(final Exception err)
			{
			return wrapException(err);
			}
		return super.initializeKnime();
		}

	
	@Override
	public void disposeKnime()
		{
		CloserUtil.close(this.bedReader);
		this.bedReader = null;
		this.intervalTreeMap=null;
		this.parsedFormat = null;
		super.disposeKnime();
		}
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		return doVcfToVcf(inputName);
		}

	
	
	public static void main(String[] args) throws Exception
		{
		new VCFBed().instanceMainWithExit(args);
		}
}
