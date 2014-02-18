package com.github.lindenb.jvarkit.util.vcf;

import java.io.File;
import java.io.IOException;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.logging.Logger;

import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.BlockCompressedOutputStream;

import org.broad.tribble.readers.LineIterator;
import org.broad.tribble.readers.LineReader;
import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFConstants;
import org.broadinstitute.variant.vcf.VCFContigHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.io.IOUtils;

public class VCFUtils
	{
	private static final Logger LOG=Logger.getLogger("jvarkit");
	
	public static class CodecAndHeader
		{
		public VCFCodec codec;
		public VCFHeader header;
		}
	/*
	private static class LR implements LineReader
		{
		LinkedList<String> stack=null;
		LR(LinkedList<String> stack)
			{
			this.stack=stack;
			}
		@Override
		public String readLine() throws IOException {
			return (stack.isEmpty()?null:stack.removeFirst());
			}
		@Override
		public void close() {
			
			}
		}*/
	private static class LIT implements LineIterator
		{
		LinkedList<String> stack=null;
		LIT(LinkedList<String> stack)
			{
			this.stack=stack;
			}
		@Override
		public boolean hasNext() {
			return !stack.isEmpty();
			}
		@Override
		public String next() {
			return stack.removeFirst();
			}
		@Override
		public String peek() {
			return stack.peek();
			}
		@Override
		public void remove() {
			throw new UnsupportedOperationException();
			}
		}
	
	public static CodecAndHeader parseHeader(LineReader r) throws IOException
		{
		CodecAndHeader vh=new CodecAndHeader();
		vh.codec=new VCFCodec();
		LinkedList<String> stack=new LinkedList<String>();
    	String line;
    	while((line=r.readLine())!=null && line.startsWith("#"))
    		{
    		stack.add(line);
    		if(line.startsWith("#CHROM\t")) break;
    		}
    	vh.header=  (VCFHeader)vh.codec.readActualHeader(new LIT(stack));
    	return vh;
		}
	
	public static CodecAndHeader parseHeader(LineIterator r)
		{
		CodecAndHeader vh=new CodecAndHeader();
		vh.codec=new VCFCodec();
		LinkedList<String> stack=new LinkedList<String>();
		while(r.hasNext())
			{
			String line=r.peek();
			if(!line.startsWith("#")) break;
			stack.add(r.next());
			if(line.startsWith("#CHROM\t")) break;
			}
		vh.header=  (VCFHeader)vh.codec.readActualHeader(new LIT(stack));
    	return vh;
		}

	public static CodecAndHeader parseHeader(List<String> list)
		{
		CodecAndHeader vh=new CodecAndHeader();
		vh.codec=new VCFCodec();
		vh.header=  (VCFHeader)vh.codec.readActualHeader(new LIT(new LinkedList<String>(list)));
		return vh;
		}
	
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
			return createVcfIteratorStdin();
			}
		else
			{
			LOG.info("reading from "+IN);
			return new VcfIterator(IOUtils.openURIForReading(IN));
			}
		}
	/** create a VCF iterator
	 * 
	 * @param IN : input uri or null for stdin
	 * */
	public static  VcfIterator createVcfIteratorStdin() throws IOException
		{
		LOG.info("reading from stdin");
		return new VcfIterator(System.in);
		}
	
	public static  VariantContextWriter createVariantContextWriterToStdout()
		{
		return VariantContextWriterFactory.create(System.out,null,EnumSet.noneOf(Options.class));
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
			return createVariantContextWriterToStdout();
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
	
	public static SAMSequenceRecord contigLineToSamSequenceRecord(String line)
		{
		if(!line.startsWith(VCFConstants.CONTIG_HEADER_START+"="))
			{
			throw new IllegalArgumentException("not a contig line "+line);
			}
		String sequencename=null;
		Integer sequencelength=null;
		String tokens1[]=line.split("[<=,>\'\"]+");
		for(int i=0;i+1< tokens1.length;++i)
			{
			if(tokens1[i].equals("ID"))
				{
				sequencename=tokens1[i+1];
				}
			else if(tokens1[i].equals("length"))
				{
				sequencelength=Integer.parseInt(tokens1[i+1]);
				}
			}
		if(sequencelength==null) throw new IllegalArgumentException("no 'length' in  contig line "+line);
		if(sequencename==null) throw new IllegalArgumentException("no 'name' in  contig line "+line);
		return new SAMSequenceRecord(sequencename, sequencelength);
		}
	
	public static String samSequenceRecordToVcfContigLine(final SAMSequenceRecord ssr)
		{
		String as=ssr.getAssembly();
		if(as==null) as="";
		as=as.trim();
		return VCFConstants.CONTIG_HEADER_START+
				"=<ID="+ssr.getSequenceName()+
				",length="+ssr.getSequenceLength()+
				(as.isEmpty()?"":",assembly="+as)+
				">";
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
			String as=ssr.getAssembly();
			if(as!=null && !as.trim().isEmpty()) mapping.put("assembly",as);
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
