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
package com.github.lindenb.jvarkit.util.vcf;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.logging.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineReader;
import htsjdk.tribble.util.TabixUtils;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCF3Codec;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;

import com.github.lindenb.jvarkit.io.IOUtils;

public class VCFUtils
	{
	private static final Logger LOG=Logger.getLogger("jvarkit");
	
	public static class CodecAndHeader
		{
		public AbstractVCFCodec codec;
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
	
	public static List<String> parseHeaderLines(LineReader r) throws IOException
		{
		LinkedList<String> stack=new LinkedList<String>();
		String line;
		while((line=r.readLine())!=null && line.startsWith("#"))
			{
			stack.add(line);
			if(line.startsWith("#CHROM\t")) break;
			}
		return stack;
		}

	
	public static CodecAndHeader parseHeader(LineReader r) throws IOException
		{
		return parseHeader(parseHeaderLines(r));
		}
	
	public static List<String> parseHeaderLines(BufferedReader r) throws IOException
		{
		LinkedList<String> stack=new LinkedList<String>();
		String line;
		while((line=r.readLine())!=null && line.startsWith("#"))
			{
			stack.add(line);
			if(line.startsWith("#CHROM\t")) break;
			}
		return stack;
		}

	
	public static CodecAndHeader parseHeader(BufferedReader r) throws IOException
		{
		return parseHeader(parseHeaderLines(r));
		}
		
	
	public static CodecAndHeader parseHeader(LineIterator r)
		{
		CodecAndHeader vh=new CodecAndHeader();
		vh.codec=null;
		LinkedList<String> stack=new LinkedList<String>();
		while(r.hasNext())
			{
			String line=r.peek();
			if(!line.startsWith("#")) break;
			stack.add(r.next());
			if(line.startsWith("#CHROM\t")) break;
			}
		vh.codec = findCodecFromLines(stack);
		vh.header=  (VCFHeader)vh.codec.readActualHeader(new LIT(stack));
    	return vh;
		}
	/** find a codec from the lines header. if not found, return default codec */
	public static AbstractVCFCodec findCodecFromLines(final List<String> list)
		{
		for(String line: list)
			{
			String formatString = line;
			if(formatString.startsWith("##"))
				{
				formatString=formatString.substring(2);
				}
			int eq= formatString.indexOf('=');
			if(eq==-1) continue;
			
			if(!VCFHeaderVersion.isFormatString(formatString.substring(0,eq))) continue;
			
			VCFHeaderVersion version=VCFHeaderVersion.getHeaderVersion(line	);
			if(version==null) continue;
			switch(version)
				{
				case VCF3_2: 
				case VCF3_3: return new VCF3Codec();
				case VCF4_0:
				case VCF4_1:
				case VCF4_2: return new VCFCodec();
				}
			}
		return createDefaultVCFCodec();
		}
	
	public static CodecAndHeader parseHeader(List<String> list)
		{
		CodecAndHeader vh=new CodecAndHeader();
		vh.codec=findCodecFromLines(list);
		vh.header=  (VCFHeader)vh.codec.readActualHeader(new LIT(new LinkedList<String>(list)));
		return vh;
		}
	
	/** convert a VCF header to line iterator. Created for serialization */
	public static LineIterator convertVCFHeaderToLineIterator(VCFHeader header)
		{	
		ByteArrayOutputStream baos=new ByteArrayOutputStream(1000);
		VariantContextWriter vcw=createVariantContextWriterToOutputStream(baos);
		vcw.writeHeader(header);
		vcw.close();
		return new LIT(new LinkedList<String>(
				Arrays.asList(
						new String(baos.toByteArray()).split("\n")
				)));
		}
	
	/** create a default VCF codec */
	public static AbstractVCFCodec createDefaultVCFCodec()
		{
		return new VCFCodec();
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
		return createVariantContextWriterToOutputStream(System.out);
		}
	
	public static  VariantContextWriter createVariantContextWriterToOutputStream(OutputStream ostream)
		{
		VariantContextWriterBuilder vcwb=new VariantContextWriterBuilder();
		vcwb.setCreateMD5(false);
		vcwb.setOutputStream(ostream);
		vcwb.setReferenceDictionary(null);
		vcwb.clearOptions();
		return vcwb.build();
		}

	
	/**
	 * create a VariantContextWriter
	 * @param OUT output file or null to stdout
	 * @return
	 * @throws IOException
	 */
	public static  VariantContextWriter createVariantContextWriter(File OUT) throws IOException
		{
		VariantContextWriterBuilder vcwb=new VariantContextWriterBuilder();
		vcwb.setCreateMD5(false);
		vcwb.setReferenceDictionary(null);
		vcwb.clearOptions();
		
		if(OUT==null)
			{
			LOG.info("writing to stdout");
			return createVariantContextWriterToStdout();
			}
		else
			{
			IOUtil.assertFileIsWritable(OUT);
			vcwb.setOutputFile(OUT);
			return vcwb.build();
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
	public static boolean isVcfFile(File f)
		{
		if(f==null || !f.isFile()) return false;
		String s=f.getName();
		return s.endsWith(".vcf") || s.endsWith(".vcf.gz");
		}
	
	
	public static boolean isTabixVcfFile(File f)
		{
		if(!isVcfFile(f)) return false;
		String filename=f.getName();
		if(!filename.endsWith(".gz")) return false;
		File index=new File(f.getParentFile(),
				filename+ TabixUtils.STANDARD_INDEX_EXTENSION
				);
		return index.exists() &&  index.isFile();
		}
	
	public static String findChromNameEquivalent(String chromName,VCFHeader h)
		{
		SAMSequenceDictionary dict=h.getSequenceDictionary();
		if(dict==null || dict.getSequence(chromName)!=null) return chromName;
		if(chromName.startsWith("chr"))
			{
			SAMSequenceRecord ssr=dict.getSequence(chromName.substring(3));
			if(ssr!=null) return ssr.getSequenceName();
			}
		else
			{
			SAMSequenceRecord ssr=dict.getSequence("chr"+chromName);
			if(ssr!=null) return ssr.getSequenceName();
			}
		
		if(chromName.equals("MT") &&  dict.getSequence("chrM")!=null) return "chrM";
		if(chromName.equals("chrM") &&  dict.getSequence("MT")!=null) return "MT";
		return null;
		}
	
	/** escape String for a VALUE in the INFO column 
	 * From the spec : no white-space, semi-colons, or equals-signs permitted; 
	 * commas are permitted only as delimiters for lists of values.
	 * INFO Fields are encoded as a semicolon-separated series of short keys with optional values in the format:
	 */
	public static String escapeInfoField(String infoValue)
		{
		return infoValue.replaceAll("[ \t\n;=,_]+", "_");
		}
	
	/** used to return any VCF Attribute as a java.util.List */
	@SuppressWarnings("unchecked")
	public static List<? extends Object> attributeAsList(Object o)
		{
		if(o==null) return Collections.emptyList();
		if(o instanceof List) return ((List<Object>)o);
		if(o.getClass().isArray())
			{	
			return Arrays.asList((Object[])o);
			}
		return Collections.singletonList(o);
		}
	
	/**
	 * Creates a comparator based on CHROM/POS/REF
	 * @return
	 */
	public static Comparator<VariantContext> createChromPosRefComparator()
		{
		return new Comparator<VariantContext>()
			{
			@Override
			public int compare(final VariantContext ctx1, final VariantContext ctx2)
				{
				int i = ctx1.getChr().compareTo(ctx2.getChr());
				if(i!=0) return i;
				int i0 = ctx1.getStart();
				int i1 = ctx2.getStart();
				if(i0 != i1) return i0-i1;
				
				i = ctx1.getReference().compareTo(
					ctx2.getReference()); 
				if(i!=0) return i;
				return 0;
				}
			};
		}
	/**
	 * Creates a comparator based on dict(CHROM)/POS/REF
	 * @return
	 */
	public static Comparator<VariantContext> createTidPosRefComparator(SAMSequenceDictionary dict)
		{
		return new _DictCompareCtx(dict);
		}

	
	private static class _DictCompareCtx
		implements Comparator<VariantContext>
		{
		SAMSequenceDictionary dict;
		_DictCompareCtx(SAMSequenceDictionary dict) {this.dict=dict;}
		
		private int tid(String chrom)
			{
			int t= dict.getSequenceIndex(chrom);
			if(t==-1) throw new IllegalArgumentException("chromosome \""+chrom+"\" is missing in dictionary");
			return t;
			}
		@Override
		public int compare(
				final VariantContext ctx1,
				final VariantContext ctx2
				)
			{			
			int i0 = tid(ctx1.getChr());
			int i1 = tid(ctx2.getChr());
			if(i0 != i1) return i0-i1;
				
			i0 = ctx1.getStart();
			i1 = ctx2.getStart();
			if(i0 != i1) return i0-i1;
			
			int i = ctx1.getReference().compareTo(ctx2.getReference()); 
			if(i!=0) return i;
			return 0;
			}

		}
	
	/** return true if 'f' is a file, path ends with '.gz' and there is an associated .tbi file */
    public static final boolean isValidTabixFile(File f)
		 	{
			if (!f.isFile())
				return false;
			String filename = f.getName();
			if (!filename.endsWith(".gz"))
				return false;
			File index = new File(f.getParentFile(), filename +TabixUtils.STANDARD_INDEX_EXTENSION);
			return index.exists() && index.isFile();
			}
    
    public static FeatureCodec<VariantContext, LineIterator> createAsciiFeatureCodec()
    	{
    	return new VariantContextCodec();
    	}
    
    /** variant feature codec for tribble index */
    private static class VariantContextCodec
		extends AsciiFeatureCodec<VariantContext>
		{
		private VCFUtils.CodecAndHeader codecHeader;
		VariantContextCodec()
			{
			super(VariantContext.class);
			}
	
		@Override
		public VariantContext decode(String line)
			{
			return this.codecHeader.codec.decode(line);
			}
		
		@Override
		public VCFHeader readActualHeader(LineIterator reader) {
			this.codecHeader=VCFUtils.parseHeader(reader);
			return this.codecHeader.header;
			}
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
