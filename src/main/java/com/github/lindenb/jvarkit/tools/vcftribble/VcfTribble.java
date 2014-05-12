package com.github.lindenb.jvarkit.tools.vcftribble;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.github.lindenb.jvarkit.util.picard.cmdline.Option;
import com.github.lindenb.jvarkit.util.picard.cmdline.Usage;

import htsjdk.samtools.util.LocationAware;
import htsjdk.samtools.util.Log;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import htsjdk.tribble.AbstractFeatureCodec;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.TribbleIndexedFeatureReader;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.AsciiLineReaderIterator;
import htsjdk.tribble.readers.LineReader;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter;




@SuppressWarnings("unchecked")
public class VcfTribble extends AbstractVCFFilter
	{
	 private static Log LOG=Log.getInstance(VcfTribble.class); 
	
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+"Get the INFO from a VCF and use it for another VCF. ";

	@Option(shortName="TB",doc="File having chrom/start/end data")
	public File TRIBBLE;
	
	@Option(shortName="FMT",doc="format. Field with ${number} will be replaced with the column of the BED.",optional=true)
	public String FORMAT="${1}:${2}-${3}";
	
	@Option(shortName="T",doc="Key for the INFO field",optional=true)
	public String TAG="TAG";
	
	private char columnDelim='\t';
	private int column_chrom=0;
	private int column_start=1;
	private int column_end=2;
	
	
	private static abstract class Chunk
		{
		public abstract String toString(List<String> tokens);
		public Chunk next=null;
		}
	
	private static class PlainChunk extends Chunk
		{
		String s;
		PlainChunk(String s){this.s=s;}
		public String toString(List<String> tokens)
			{
			return s+(next==null?"":next.toString(tokens));
			}
		}
	private static class ColChunk extends Chunk
		{
		int index;
		ColChunk(int index){ this.index=index;}
		public String toString(List<String> tokens)
			{
			return tokens.get(index)+(next==null?"":next.toString(tokens));
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

	
	private class MyFeatureCodec
		extends AbstractFeatureCodec<MyFeature, PositionalBufferedStream>
		{
		protected MyFeatureCodec() {
			super(MyFeature.class);
			}
		@Override
		public boolean canDecode(String filename)
			{
			return true;
			}
		private List<String> nextLine(PositionalBufferedStream in)
			throws IOException
			{
			List<String> tokens=null;
			StringBuilder currStr=null;
			int c;
			while((c=in.read())!=-1 && c!='\n')
				{
				if(c==columnDelim)
					{
					if(tokens==null)
						{
						tokens=new ArrayList<String>();
						tokens.add("");//first character was a delimiter
						}
					if(currStr!=null) tokens.add(currStr.toString());
					currStr=null;
					continue;
					}
				if(currStr!=null) currStr=new StringBuilder();
				currStr.append((char)c);
				}
			if(currStr!=null)
				{
				if(tokens==null) tokens=new ArrayList<String>();
				tokens.add(currStr.toString());
				}
			return tokens;
			}
		
		
		
		@Override
		public MyFeature decode(PositionalBufferedStream in)
				throws IOException
			{
			List<String> tokens;
			while((tokens=nextLine(in))!=null)
				{
				if(tokens.isEmpty()) continue;
				if(tokens.get(0).startsWith("#")) continue;
				return new MyFeature(tokens);
				}
			return null;
			}
		
		public Feature decodeLoc(PositionalBufferedStream in) throws IOException
			{
			return decode(in);
			}
		@Override
		public Class<MyFeature> getFeatureType()
			{
			return MyFeature.class;
			}
		@Override
		public FeatureCodecHeader readHeader(PositionalBufferedStream in)
				throws IOException
			{
			StringBuilder b=new StringBuilder();
			for(;;)
				{
				int c=in.peek();
				if(c==-1  || c!='#') break;
				while((c=in.read())!=-1 && c!='\n')
					{
					b.append((char)c);
					}
				}
			return new FeatureCodecHeader(b.toString(),in.getPosition());
			}
		@Override
		public PositionalBufferedStream makeSourceFromStream(
				InputStream bufferedInputStream) {

	        final PositionalBufferedStream pbs;
	        if (bufferedInputStream instanceof PositionalBufferedStream) {
	           return (PositionalBufferedStream) bufferedInputStream;
	        } else {
	            return new PositionalBufferedStream(bufferedInputStream);
	        }
			}
		@Override
		public LocationAware makeIndexableSourceFromStream(
				InputStream bufferedInputStream) {
	        final PositionalBufferedStream pbs;
	        if (bufferedInputStream instanceof PositionalBufferedStream) {
	           return (PositionalBufferedStream) bufferedInputStream;
	        } else {
	            return new PositionalBufferedStream(bufferedInputStream);
	        }
			}
		@Override
		public boolean isDone(PositionalBufferedStream source) {
			try {
				return source.isDone();
			} catch (IOException e) {
				return true;
			}
		}
		@Override
		public void close(PositionalBufferedStream source) {
			source.close();
			
		}
		}
	
	private class MyFeature
		implements Feature
		{
		List<String> tokens;
		
		MyFeature(List<String> tokens)
			{
			this.tokens=tokens;
			}
		
		@Override
		public String getChr()
			{
			return this.tokens.get(column_chrom);
			}
		
		@Override
		public int getStart()
			{
			return Integer.parseInt(this.tokens.get(column_start));
			}
		
		@Override
		public int getEnd()
			{
			return  Integer.parseInt(this.tokens.get(column_end));
			}
		
		}
	
	@SuppressWarnings("rawtypes")
	@Override
	protected void doWork(VcfIterator r, VariantContextWriter w)
			throws IOException
		{
		VCFCodec codeIn1=new VCFCodec();	
		String line;
		Chunk parsedFormat=parseFormat(this.FORMAT);
		if(parsedFormat==null)parsedFormat=new PlainChunk("");

		LOG.info("opening "+this.TRIBBLE);
		
		MyFeatureCodec tribbleCodec=new MyFeatureCodec();
		
		 Index tribbleIndex=null;
		 File idxFile = Tribble.indexFile(this.TRIBBLE);
		 if(idxFile.exists() && idxFile.lastModified()>=this.TRIBBLE.lastModified())
		 	{
			LOG.info("loading index in memory for "+this.TRIBBLE);
			tribbleIndex=IndexFactory.createLinearIndex(TRIBBLE, tribbleCodec);
		 	}
		 else
		 	{
			LOG.info("loading index from file "+this.TRIBBLE);
			tribbleIndex=IndexFactory.loadIndex(idxFile.getPath());
		 	}
		
		 
		
	    TribbleIndexedFeatureReader<MyFeature,PositionalBufferedStream> tribbleReader= new TribbleIndexedFeatureReader(
	    		this.TRIBBLE.getPath(),
	    		tribbleCodec,
	    		tribbleIndex
	    		);
	   
		VCFHeader header1=r.getHeader();
		
		VCFHeader h2=new VCFHeader(header1.getMetaDataInInputOrder(),header1.getSampleNamesInOrder());
		h2.addMetaDataLine(new VCFInfoHeaderLine(TAG, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "metadata added from "+TRIBBLE+" . Format was "+FORMAT));
		
		w.writeHeader(h2);
		while(r.hasNext())
			{
			VariantContext ctx=r.next();
			Set<String> annotations=new HashSet<String>();

			
			
			
			CloseableTribbleIterator<MyFeature> iter= tribbleReader.query(
					ctx.getChr(),
					ctx.getStart(),
					ctx.getEnd()
					);
			while(iter.hasNext())
				{
				MyFeature feat=iter.next();
				String newannot=parsedFormat.toString(feat.tokens);
				if(!newannot.isEmpty())
					annotations.add(newannot.replaceAll("[ , ;=]","_"));
				
				if(annotations.isEmpty())
					{
					w.add(ctx);
					continue;
					}
				}
			iter.close();
			
		
			VariantContextBuilder vcb=new VariantContextBuilder(ctx);
			vcb.attribute(TAG, annotations.toArray());
			w.add(vcb.make());
			}
		tribbleReader.close();
		}
	
	
	public static void main(String[] args) throws IOException
		{
		new VcfTribble().instanceMainWithExit(args);
		}
}
