package com.github.lindenb.jvarkit.tools.vcftribble;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;

import org.broad.tribble.CloseableTribbleIterator;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.FeatureCodecHeader;
import org.broad.tribble.TribbleIndexedFeatureReader;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.LineReader;
import org.broad.tribble.readers.PositionalBufferedStream;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.tools.vcfbed.VCFBed.Chunk;
import com.github.lindenb.jvarkit.tools.vcfbed.VCFBed.ColChunk;
import com.github.lindenb.jvarkit.tools.vcfbed.VCFBed.PlainChunk;
import com.github.lindenb.jvarkit.util.AbstractVCFFilter;





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

	
	private class MyFeatureCodec
		implements FeatureCodec<MyFeature>
		{
		@Override
		public boolean canDecode(String arg0)
			{
			// TODO Auto-generated method stub
			return false;
			}
		@Override
		public MyFeature decode(PositionalBufferedStream in)
				throws IOException
			{
			// TODO Auto-generated method stub
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
		public FeatureCodecHeader readHeader(PositionalBufferedStream arg0)
				throws IOException
			{
			// TODO Auto-generated method stub
			return null;
			}
		}
	
	private class MyFeature
		implements Feature
		{
		String tokens[];
		@Override
		public String getChr()
			{
			// TODO Auto-generated method stub
			return null;
			}
		
		@Override
		public int getStart()
			{
			// TODO Auto-generated method stub
			return 0;
			}
		
		@Override
		public int getEnd()
			{
			// TODO Auto-generated method stub
			return 0;
			}
		
		}
	
	@Override
	protected void doWork(LineReader in, VariantContextWriter w)
			throws IOException
		{
		VCFCodec codeIn1=new VCFCodec();	
		VCFCodec codeIn3=new VCFCodec();	
		String line;
		Chunk parsedFormat=parseFormat(this.FORMAT);
		if(parsedFormat==null)parsedFormat=new PlainChunk("");

		StringWriter sw=new StringWriter();
		LOG.info("opening "+this.TRIBBLE);
		
		MyFeatureCodec tribbleCodec=new MyFeatureCodec();
		
		Index tribbleIndex=IndexFactory.createLinearIndex(TRIBBLE, tribbleCodec);
		
	    TribbleIndexedFeatureReader<MyFeature> tribbleReader= new TribbleIndexedFeatureReader(this.TRIBBLE.getPath(),tribbleCodec,tribbleIndex);
	   
		VCFHeader header1=(VCFHeader)codeIn1.readHeader(in);
		
		VCFHeader h2=new VCFHeader(header1.getMetaDataInInputOrder(),header1.getSampleNamesInOrder());

		
		w.writeHeader(h2);
		while((line=in.readLine())!=null)
			{
			VariantContext ctx=codeIn1.decode(line);
			Set<String> annotations=new HashSet<String>();

			
			
			List<VariantContext> variantsList=new ArrayList<VariantContext>();
			CloseableTribbleIterator<MyFeature> iter= (CloseableTribbleIterator<MyFeature>)tribbleReader.query(
					ctx.getChr(),
					ctx.getStart(),
					ctx.getEnd()
					);
			while(iter.hasNext())
				{
				MyFeature feat=iter.next();
				String tokens[]=feat.tokens;
				String newannot=parsedFormat.toString(tokens);
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
