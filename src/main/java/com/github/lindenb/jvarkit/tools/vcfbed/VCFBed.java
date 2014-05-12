package com.github.lindenb.jvarkit.tools.vcfbed;

import java.io.IOException;
import java.util.HashSet;
import java.util.Set;
import java.util.regex.Pattern;


import htsjdk.samtools.cmdline.Option;
import htsjdk.samtools.cmdline.Usage;
import htsjdk.samtools.util.Log;


import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


/**
 * 
 * VCFBed
 *
 */
public class VCFBed extends AbstractVCFFilter
	{
	 private static Log LOG=Log.getInstance(VCFBed.class);
	
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+" Cross information between a VCF and a BED .";

	@Option(shortName="FMT",doc="format. Field with ${number} will be replaced with the column of the BED.",optional=true)
	public String FORMAT="${1}:${2}-${3}";
	
	@Option(shortName="TBX",doc="BED file indexed with tabix",optional=false)
	public String TABIX;
	
	@Option(shortName="T",doc="Key for the INFO field",optional=true)
	public String TAG="TAG";
	
	
	private static abstract class Chunk
		{
		public abstract String toString(String tokens[]);
		public Chunk next=null;
		}
	
	private static class PlainChunk extends Chunk
		{
		String s;
		PlainChunk(String s){this.s=s;}
		public String toString(String tokens[])
			{
			return s+(next==null?"":next.toString(tokens));
			}
		}
	private static class ColChunk extends Chunk
		{
		int index;
		ColChunk(int index){ this.index=index;}
		public String toString(String tokens[])
			{
			return tokens[index]+(next==null?"":next.toString(tokens));
			}
		}

	
	private Chunk parseFormat(String s)
		{
		if(s==null || s.isEmpty()) return null;
		if(s.startsWith("${"))
			{
			int j=s.indexOf('}',2);
			if(j==-1) throw new IllegalArgumentException("bad format in \""+s+"\".");
			try
				{
				int col=Integer.parseInt(s.substring(2, j).trim());
				if(col<1) throw new IllegalArgumentException();
				ColChunk c=new ColChunk(col-1);
				c.next=parseFormat(s.substring(j+1));
				return c;
				}
			catch(Exception err)
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
		StringBuilder sb=new StringBuilder();
		while(i< s.length() && s.charAt(i)!='$')
			{
			sb.append(s.charAt(i));
			i++;
			}
		PlainChunk c=new PlainChunk(sb.toString());
		c.next=parseFormat(s.substring(i));
		return c;
		}
	
	@Override
	public String getVersion()
		{
		return "1.0";
		}
	
	@Override
	protected void doWork(VcfIterator r, VariantContextWriter w)
			throws IOException
		{
		Pattern tab=Pattern.compile("[\t]");
		LOG.info("parsing "+this.FORMAT);
		Chunk parsedFormat=parseFormat(this.FORMAT);
		if(parsedFormat==null)parsedFormat=new PlainChunk("");
		
		LOG.info("opening TABIX "+this.TABIX);
		TabixReader tabix= new TabixReader(this.TABIX);
		
		VCFHeader header=r.getHeader();

		VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
		h2.addMetaDataLine(new VCFInfoHeaderLine(TAG, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "metadata added from "+TABIX+" . Format was "+FORMAT));
		w.writeHeader(h2);
		while(r.hasNext())
			{
			VariantContext ctx=r.next();
			Set<String> annotations=new HashSet<String>();
			String line2;
			
			
			
			String q=ctx.getChr()+":"+ctx.getStart()+"-"+(ctx.getEnd());
			int reg[]=tabix.parseReg(q);
			if(reg[0]==-1)
				{
				w.add(ctx);
				continue;
				}
			TabixReader.Iterator iter=tabix.query(reg[0],reg[1],reg[2]);
			while(iter!=null && (line2=iter.next())!=null)
				{
				String tokens[]=tab.split(line2);
				String newannot=parsedFormat.toString(tokens);
				if(!newannot.isEmpty())
					annotations.add(newannot.replaceAll("[ , ;=]+","_"));
				}
			if(annotations.isEmpty())
				{
				w.add(ctx);
				continue;
				}
			VariantContextBuilder vcb=new VariantContextBuilder(ctx);
			vcb.attribute(TAG, annotations.toArray());
			w.add(vcb.make());
			}
		tabix.close();
		}
	
	public static void main(String[] args) throws Exception
		{
		new VCFBed().instanceMainWithExit(args);
		}
}
