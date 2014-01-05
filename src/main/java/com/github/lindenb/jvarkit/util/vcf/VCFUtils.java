package com.github.lindenb.jvarkit.util.vcf;

import java.io.File;
import java.io.IOException;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.logging.Logger;

import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.BlockCompressedOutputStream;

import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;
import org.broadinstitute.variant.vcf.VCFContigHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;

public class VCFUtils
	{
	private static final Logger LOG=Logger.getLogger("jvarkit");
	
	
	/** create a VCF iterator
	 * 
	 * @param IN input vcf file
	 * */
	public static  VcfIterator createVcfIteratorFromFile(File f) throws IOException
		{
		LOG.info("reading vcf from "+f);
		return new VcfIterator(IOUtils.openFileForReading(f));	
		}
	
	/** create a VCF iterator
	 * 
	 * @param IN : input uri or null for stdin
	 * */
	public static  VcfIterator createVcfIterator(String IN) throws IOException
		{
		if(IN==null)
			{
			LOG.info("reading from stdin");
			return new VcfIterator(System.in);
			}
		else
			{
			LOG.info("reading from "+IN);
			return new VcfIterator(IOUtils.openURIForReading(IN));
			}
		}

	/**
	 * create a VariantContextWriter
	 * @param OUT output file or null to stdout
	 * @return
	 * @throws IOException
	 */
	public static  VariantContextWriter createVariantContextWriter(File OUT) throws IOException
		{
		if(OUT==null)
			{
			LOG.info("writing to stdout");
			return VariantContextWriterFactory.create(System.out,null,EnumSet.noneOf(Options.class));
			}
		else if(OUT.getName().endsWith(".gz"))
			{
			LOG.info("writing to "+OUT+" as bgz file.");
			BlockCompressedOutputStream bcos=new BlockCompressedOutputStream(OUT);
			return VariantContextWriterFactory.create(bcos,null,EnumSet.noneOf(Options.class));
			}
		else
			{
			LOG.info("writing to "+OUT);
			return  VariantContextWriterFactory.create(OUT,null,EnumSet.noneOf(Options.class));
			}
		}
	
	public static SortedSet<VCFContigHeaderLine>
		samSequenceDictToVCFContigHeaderLine(SAMSequenceDictionary dict)
		{
		SortedSet<VCFContigHeaderLine> meta2=new TreeSet<VCFContigHeaderLine>();
		for(SAMSequenceRecord ssr: dict.getSequences())
			{
			Map<String,String> mapping=new HashMap<String,String>();
			mapping.put("ID", ssr.getSequenceName());
			mapping.put("length",String.valueOf(ssr.getSequenceLength()));
			VCFContigHeaderLine h=new VCFContigHeaderLine(mapping,ssr.getSequenceIndex());
			meta2.add(h);
			}
		return meta2;
		}
	
	/*
	private Pattern semicolon=Pattern.compile("[;]");
	public VCFUtils()
		{
		
		}
	
	
	public String joinInfo( Map<String,String> infoMap)
		{
		StringBuilder b=new StringBuilder();
		for(String key:infoMap.keySet())
			{
			if(isEmpty(key)) continue;
			String v=infoMap.get(key);
			if(b.length()!=0) b.append(";");
			if(isEmpty(v))
				{
				b.append(key);
				}
			else
				{
				b.append(key).append("=").append(v);
				}
			}
		return b.toString();
		}
	
	public String joinFilters( Set<String> s)
		{
		StringBuilder b=new StringBuilder();
		for(String f:s)
			{
			if(isEmpty(f)) continue;
			if(b.length()!=0) b.append(";");
			b.append(f);
			}
		if(b.length()==0) return ".";
		return b.toString();
		}
	
	public Set<String> parseFilters(String filterField)
		{
		Set<String> S=new LinkedHashSet<String>();
		if(isEmpty(filterField)) return S;
		for(String f:semicolon.split(filterField))
			{
			if(isEmpty(f)) continue;
			S.add(f);
			}
		return S;
		}
		
	
	public Map<String,String> parseInfo(String infoField)
		{
		Map<String,String> m=new LinkedHashMap<String,String>();
		if(isEmpty(infoField)) return m;
		for(String info:semicolon.split(infoField))
			{
			if(info.isEmpty()) continue;
			int eq=info.indexOf('=');
			String key;
			String value;
			if(eq!=-1)
				{
				key=info.substring(0,eq);
				value=info.substring(eq+1);
				}
			else
				{
				key=info;
				value="";
				}
			m.put(key,value);
			}
	
		return m;
		}
	
	public boolean isEmpty(String s)
		{
		return s==null || s.equals(".") || s.isEmpty();
		}
		*/
	}
