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
* 2014 creation

*/
package com.github.lindenb.jvarkit.util.vcf;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.SortedSet;
import java.util.TreeSet;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineReader;
import htsjdk.tribble.util.TabixUtils;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCF3Codec;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderVersion;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;


import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

public class VCFUtils
	{	
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
	
	public static List<String> parseHeaderLines(final LineReader r) throws IOException
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

	
	public static CodecAndHeader parseHeader(final LineReader r) throws IOException
		{
		return parseHeader(parseHeaderLines(r));
		}
	
	public static List<String> parseHeaderLines(final BufferedReader r) throws IOException
		{
		final LinkedList<String> stack=new LinkedList<String>();
		String line;
		while((line=r.readLine())!=null && line.startsWith("#"))
			{
			stack.add(line);
			if(line.startsWith("#CHROM\t")) break;
			}
		return stack;
		}

	
	public static CodecAndHeader parseHeader(final BufferedReader r) throws IOException
		{
		return parseHeader(parseHeaderLines(r));
		}
		
	
	public static CodecAndHeader parseHeader(final LineIterator r)
		{
		final CodecAndHeader vh=new CodecAndHeader();
		vh.codec=null;
		LinkedList<String> stack=new LinkedList<String>();
		while(r.hasNext())
			{
			final String line=r.peek();
			if(!line.startsWith("#")) break;
			stack.add(r.next());
			if(line.startsWith("#CHROM\t")) break;
			}
		vh.codec = findCodecFromLines(stack);
		vh.header=  (VCFHeader)vh.codec.readActualHeader(new LIT(stack));
    	return vh;
		}
	
	
	/** stringent insertion in VCF header */
	public static void safeAddMetaDataHeaderLine(final VCFHeader header,final VCFHeaderLine hl){
		final VCFHeaderLine prev =header.getMetaDataLine(hl.getKey());
		if(prev!=null && !prev.getValue().equals(hl.getValue())) {
			throw new IllegalArgumentException("Cannot insert vcf header "+hl.getKey()+" because it is already defined in VCF header");
		}
		header.addMetaDataLine(hl);
	}

	/** find a codec from the lines header. if not found, return default codec */
	public static AbstractVCFCodec findCodecFromLines(final List<String> list)
		{
		for(final String line: list)
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
	
	public static CodecAndHeader parseHeader(final List<String> list)
		{
		final CodecAndHeader vh=new CodecAndHeader();
		vh.codec=findCodecFromLines(list);
		vh.header=  (VCFHeader)vh.codec.readActualHeader(new LIT(new LinkedList<String>(list)));
		return vh;
		}
	
	/** convert a VCF header to List of String. Created for serialization */
	public static List<String> convertVCFHeaderToList(final VCFHeader header)
		{	
		final ByteArrayOutputStream baos=new ByteArrayOutputStream(1000);
		final VariantContextWriter vcw=createVariantContextWriterToOutputStream(baos);
		vcw.writeHeader(header);
		vcw.close();
		return  Arrays.asList(
						new String(baos.toByteArray()).split("\n")
				);
		}

	
	/** convert a VCF header to line iterator. Created for serialization */
	public static LineIterator convertVCFHeaderToLineIterator(final VCFHeader header)
		{	
		return new LIT(new LinkedList<String>(convertVCFHeaderToList(header)));
		}
	
	/** create a default VCF codec */
	public static AbstractVCFCodec createDefaultVCFCodec()
		{
		return new VCFCodec();
		}
	
	/** create a VCF iterator
	 * 
	 * @param IN input stream
	 * */
	@Deprecated
	public static  VCFIterator createVCFIteratorFromStream(final InputStream in) throws IOException
		{
		return createVCFIteratorFromInputStream(in);
		}
	
	/** create a VCF iterator
	 * 
	 * @param IN input vcf file
	 * */
	public static  VCFIterator createVCFIteratorFromFile(final File vcfOrBcfFile) throws IOException
		{
		return createVCFIteratorFromPath(vcfOrBcfFile==null?null:vcfOrBcfFile.toPath());
		}
	
	
	public static  VCFIterator createVCFIteratorFromPath(final Path vcfOrBcfFile) throws IOException
		{
		IOUtil.assertFileIsReadable(vcfOrBcfFile);
		return new VCFIteratorBuilder().open(vcfOrBcfFile);
		}

	
	/** create a VCF iterator
	 * 
	 * @param IN input vcf file
	 * */
	public static  VCFIterator createVCFIteratorFromInputStream(final InputStream in) throws IOException
		{
		return new VCFIteratorBuilder().open(in);	
		}

	/** create a VCF iterator from LineReader
	 * 
	 * @param IN input vcf file
	 * */
	public static  VCFIterator createVCFIteratorFromLineIterator(
			final LineIterator lineIterator,
			boolean allowConcatenatedVcf
			) throws IOException
		{
		return new VCFIteratorLineIterator(lineIterator,allowConcatenatedVcf);	
		}

	
	/** create a VCF iterator
	 * 
	 * @param IN : input uri or null for stdin
	 * */
	public static  VCFIterator createVCFIterator(final String IN) throws IOException
		{
		if(IN==null)
			{
			return createVCFIteratorStdin();
			}
		else if(Arrays.asList(IOUtil.VCF_EXTENSIONS).stream().anyMatch(S->IN.endsWith(S))
				&& !IOUtils.isRemoteURI(IN))
			{
			final File bcfFile = new File(
					IN.startsWith("file://")?
					IN.substring(7):
					IN
					);
			return createVCFIteratorFromFile(bcfFile);
			}
		else
			{
			return createVCFIteratorFromInputStream(IOUtils.openURIForReading(IN));
			}
		}
	/** create a VCF iterator
	 * 
	 * @param IN : input uri or null for stdin
	 * */
	public static  VCFIterator createVCFIteratorStdin() throws IOException
		{
		return createVCFIteratorFromInputStream(System.in);
		}
	
	public static  VariantContextWriter createVariantContextWriterToStdout()
		{
		return createVariantContextWriterToOutputStream(System.out);
		}
	
	public static  VariantContextWriter createVariantContextWriterToOutputStream(final OutputStream ostream)
		{
		final VariantContextWriterBuilder vcwb=new VariantContextWriterBuilder();
		vcwb.setCreateMD5(false);
		vcwb.setOutputStream(ostream);
		vcwb.setReferenceDictionary(null);
		vcwb.clearOptions();
		return new VariantContextWriterDelayedFlush(vcwb.build());
		}

	/** wrap delegate into a VCF iterator printing progress */
	public static VCFIterator wrapInProgress(final VCFIterator delegate)
		{
		return new VCFIteratorProgress(delegate);
		}
	
	/** checkError but don't flush everything each time */
	private static class VariantContextWriterDelayedFlush
		extends DelegateVariantContextWriter
		{
		boolean lastCheckError=false;
		long lastCheck=System.currentTimeMillis();
		VariantContextWriterDelayedFlush(final VariantContextWriter delegate) {
			super(delegate);
		}
		
		@Override
		public boolean checkError() {
			if(!this.lastCheckError ) {
				final long now = System.currentTimeMillis();
				if((now-lastCheck)>10*1000) //10sec
					{
					this.lastCheckError =  getDelegate().checkError();
					this.lastCheck=now;
					}
				}
			return this.lastCheckError;
			}
		}
	
	/**
	 * create a VariantContextWriter
	 * @param OUT output path or null to stdout
	 * @return
	 * @throws IOException
	 */
	public static  VariantContextWriter createVariantContextWriterToPath(final Path pathorNull) throws IOException
		{
		if(pathorNull==null) {
			return createVariantContextWriterToStdout();
			}
		final VariantContextWriterBuilder vcwb=new VariantContextWriterBuilder();
		vcwb.setCreateMD5(false);
		vcwb.setReferenceDictionary(null);
		vcwb.clearOptions();
		vcwb.setOutputPath(pathorNull);
		return new VariantContextWriterDelayedFlush(vcwb.build());
		}
	
	/**
	 * create a VariantContextWriter
	 * @param OUT output file or null to stdout
	 * @return
	 * @throws IOException
	 */
	public static  VariantContextWriter createVariantContextWriter(final File OUT) throws IOException
		{
		return createVariantContextWriterToPath(OUT==null?null:OUT.toPath());
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
	public static boolean isVcfFile(final File f)
		{
		return isVcfPath(f==null?null:f.toPath());
		}
	
	/** return true if p ends with vcf/vcf.gz/vcf.bgz */
	public static boolean isVcfPath(final Path f)
		{
		if(f==null || Files.isDirectory(f)) return false;
		final String s=f.getFileName().toString();
		return StringUtils.endsWith(s,
				IOUtil.VCF_FILE_EXTENSION,
				IOUtil.COMPRESSED_VCF_FILE_EXTENSION,
				IOUtil.BCF_FILE_EXTENSION
				);
		}
	
	/** returns true if file ends with .vcf.gz and a .tbi file is associated */
	public static boolean isTabixVcfFile(final File f)
		{
		return isTabixVcfPath(f==null?null:f.toPath());
		}
	
	/** returns true if file ends with .vcf and a .idx file is associated */
	public static boolean isTribbleVcfFile(final File f)
		{
		return isTribbleVcfPath(f==null?null:f.toPath());
		}
	
	/** returns true if file ends with .vcf.gz and a .tbi file is associated */
	public static boolean isTabixVcfPath(final Path f)
		{
		if(!isVcfPath(f)) return false;
		final String filename=f.getFileName().toString();
		if(!(filename.endsWith(".vcf.gz")|| filename.endsWith(".vcf.bgz"))) return false;
		final Path index = Paths.get(
				f.getParent().toString(),
				filename+ TabixUtils.STANDARD_INDEX_EXTENSION
				);
		return Files.exists(index) && !Files.isDirectory(index);
		}
	
	/** returns true if file ends with .vcf and a .idx file is associated */
	public static boolean isTribbleVcfPath(final Path f)
		{
		if(!isVcfPath(f)) return false;
		final String filename=f.getFileName().toString();
		if(!filename.endsWith(".vcf")) return false;
		final Path index=Paths.get(
				f.getParent().toString(),
				filename+ Tribble.STANDARD_INDEX_EXTENSION
				);
		return Files.exists(index) && !Files.isDirectory(index);
		}

	
	@Deprecated//use contigNameConverter
	public static String findChromNameEquivalent(String chromName,VCFHeader h)
		{
		final SAMSequenceDictionary dict=h.getSequenceDictionary();
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
	
	/**
	 * 
	 */
	public CloseableIterator<VCFIterator> readConcatenatedVcfFile(final BufferedReader r) throws IOException {
		return new VcfFileIterator(IOUtils.toLineIterator(r));
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
				int i = ctx1.getContig().compareTo(ctx2.getContig());
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
		final ContigDictComparator contigDictComparator;
		_DictCompareCtx(final SAMSequenceDictionary dict)
			{
			this.contigDictComparator= new ContigDictComparator(dict);
			}
		
		
		@Override
		public int compare(
				final VariantContext ctx1,
				final VariantContext ctx2
				)
			{			
			int i = this.contigDictComparator.compare(ctx1.getContig(), ctx2.getContig());
			if(i!=0) return i;
				
			i = Integer.compare(ctx1.getStart() , ctx2.getStart());
			if(i!=0) return i;
			
			i = ctx1.getReference().compareTo(ctx2.getReference()); 
			if(i!=0) return i;
			return 0;
			}

		}
	
	/**
	 * Creates a comparator based on dict(CHROM)/POS/REF
	 * @return
	 */
	public static Comparator<VariantContext> createTidPosComparator(final SAMSequenceDictionary dict)
		{
		return new _DictCompareCtxTidPos(dict);
		}

	
	private static class _DictCompareCtxTidPos
		implements Comparator<VariantContext>
		{
		final SAMSequenceDictionary dict;
		_DictCompareCtxTidPos(final SAMSequenceDictionary dict) {
			this.dict=dict;
			if(this.dict==null) {
				throw new JvarkitException.VcfDictionaryMissing("No Sequence Dictionary provided.");
				}
			}
		
		private int tid(final String chrom)
			{
			int t= this.dict.getSequenceIndex(chrom);
			if(t==-1) throw new JvarkitException.ContigNotFoundInDictionary(chrom,this.dict);
			return t;
			}
		@Override
		public int compare(
				final VariantContext ctx1,
				final VariantContext ctx2
				)
			{			
			int i0 = tid(ctx1.getContig());
			int i1 = tid(ctx2.getContig());
			if(i0 != i1) return i0-i1;
				
			i0 = ctx1.getStart();
			i1 = ctx2.getStart();
			if(i0 != i1) return i0-i1;
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
	    @Override
	    public boolean canDecode(final String path) {
	        return path.toLowerCase().endsWith(".vcf");
	    }

		}

    
    public static VCFInfoHeaderLine renameVCFInfoHeaderLine(VCFInfoHeaderLine h,String name)
    	{
    	if(h.getCountType()==VCFHeaderLineCount.INTEGER)
			{
			return new VCFInfoHeaderLine(
					name,
					h.getCount(),
					h.getType(),
					h.getDescription()
					);
			}
		else
			{
			return new VCFInfoHeaderLine(
					name,
					h.getCountType(),
					h.getType(),
					h.getDescription()
					);
			}
    	}
	
    public static VCFFormatHeaderLine renameVCFFormatHeaderLine(VCFFormatHeaderLine h,String name)
		{
		if(h.getCountType()==VCFHeaderLineCount.INTEGER)
			{
			return new VCFFormatHeaderLine(
					name,
					h.getCount(),
					h.getType(),
					h.getDescription()
					);
			}
		else
			{
			return new VCFFormatHeaderLine(
					name,
					h.getCountType(),
					h.getType(),
					h.getDescription()
					);
			}
		}

    public static VCFFilterHeaderLine  renameVCFFilterHeaderLine(VCFFilterHeaderLine h,String name)
		{
    	return new VCFFilterHeaderLine(
				name,
				h.getValue()
				);
		}
    
    /** wrap a VCF iterator that checks that variants are sorted in a VCF file */
    public static VCFIterator createAssertSortedVCFIterator(final VCFIterator iter, final Comparator<VariantContext> cmp) {
    	return new AssertSortedVCFIterator(iter,cmp);
    }
    
    private static class AssertSortedVCFIterator implements VCFIterator {
    	private final VCFIterator iter;
    	private final Comparator<VariantContext> cmp;
    	private VariantContext prev = null;
    	AssertSortedVCFIterator(final VCFIterator iter, final Comparator<VariantContext> cmp) {
    		this.iter = iter;
    		this.cmp = cmp;
    		this.prev = null;
    		}
    	@Override
    	public boolean hasNext() {
    		return this.iter.hasNext();
    		}
    	@Override
    	public VariantContext peek() {
    		final VariantContext ctx =  this.iter.peek();
    		if(ctx!=null && prev!=null && this.cmp.compare(this.prev,ctx)>0) {
    			throw new RuntimeIOException("VCF is not sorted got\n  "+
    					ctx+"\nafter\n "+this.prev);
    		}
    		return ctx;
    		}
    	@Override
    	public VariantContext next() {
    		final VariantContext ctx =  this.iter.next();
    		if(ctx!=null && prev!=null && this.cmp.compare(this.prev,ctx)>0) {
    			throw new RuntimeIOException("VCF is not sorted got\n  "+
    					ctx+"\nafter\n "+this.prev);
    		}
    		this.prev = ctx;
    		return ctx;
    		}
    	//@Override
    	// public AbstractVCFCodec getCodec() { return this.iter.getCodec(); }
    	@Override
    	public VCFHeader getHeader() {
    		return this.iter.getHeader();
    		}
    	@Override
    	public void close() {
    		CloserUtil.close(this.iter);
    		}
    	}
    
    
    private static class VcfFileIterator extends AbstractIterator<VCFIterator>
    	implements CloseableIterator<VCFIterator>{
    	private  LineIterator lr;
    	public VcfFileIterator(final LineIterator lr) {
		this.lr = lr;
    	}
    	
    	@Override
    	protected VCFIterator advance() {
    		if(lr==null) return null;
    		if(!lr.hasNext()) { CloserUtil.close(lr);lr=null; return null;}
    		return new MyVCFIterator() ;
    		}
    	
    	@Override
    	public void close() {
    		if(lr==null) return;
    		CloserUtil.close(lr);
    		lr=null;
    		}
    	
    	private class MyVCFIterator extends VcfIteratorImpl
    		{
    		MyVCFIterator() {
    			super(VcfFileIterator.this.lr);
    			}
    		
    		@Override
    		public void close() {
    			if(lr==null) return;
    			while(lr.hasNext() && !lr.peek().startsWith("#"))
    				{
    				lr.next();
    				}
    			}
    		
    		}
    }
    
    private static class VCFIteratorLineIterator implements VCFIterator {
    	private final CodecAndHeader cah;
    	final LineIterator lineIter;
    	final boolean allowConcatenatedVcf;
    	private boolean closed=false;
    	VCFIteratorLineIterator(final LineIterator lineIter,
    			boolean allowConcatenatedVcf
    			)
    		{
    		this.lineIter = lineIter;
    		this.allowConcatenatedVcf = allowConcatenatedVcf;
			/* parse VCF header lines */
			final java.util.List<String> headerLines = new java.util.ArrayList<String>();
			while(lineIter.hasNext() && lineIter.peek().startsWith(VCFHeader.HEADER_INDICATOR)) {
				headerLines.add(lineIter.next());
			}
			/* parse VCF header */
			this.cah = com.github.lindenb.jvarkit.util.vcf.VCFUtils.parseHeader(headerLines);
    		}
    	
    	//@Override
    	public AbstractVCFCodec getCodec() {
    		return this.cah.codec;
    		}
    	@Override
    	public VCFHeader getHeader() {
    		return this.cah.header;
    		}
    	
    	@Override
    	public boolean hasNext() {
    		if(closed || !this.lineIter.hasNext()) return false;
    		if(this.allowConcatenatedVcf && this.lineIter.peek().startsWith(VCFHeader.HEADER_INDICATOR)) {
    			return false;
    		}
    		return true;
    		}
    	
    	@Override
    	public VariantContext next() {
    		if(!hasNext()) throw new IllegalStateException("no such variant");
    		return getCodec().decode(this.lineIter.next());
    		}
    	
    	@Override
    	public VariantContext peek() {
    		if(!hasNext()) throw new IllegalStateException("no such variant");
    		return getCodec().decode(this.lineIter.peek());
    		}
    	
    	@Override
    	public void close()  {
    		if(closed) return;
    		closed=true;
    		/* consumme remaining variants */
    		if(this.allowConcatenatedVcf) {
    			while(this.lineIter.hasNext() && !this.lineIter.peek().startsWith(VCFHeader.HEADER_INDICATOR)) {
    				this.lineIter.next();
    			}
    		} else 
    			{
    			CloserUtil.close(this.lineIter);
    			}
    		}
    	
    	
    	}
    
	@Deprecated
    private static class VCFIteratorProgress  implements VCFIterator
    	{
    	private final VCFIterator delegate;
    	private SAMSequenceDictionaryProgress progress;
    	
    	VCFIteratorProgress(final VCFIterator delegate) {
    		this.delegate=delegate;
    		this.progress = new SAMSequenceDictionaryProgress(delegate.getHeader());
    	}
    	
		@Override
		public boolean hasNext() {
			return delegate.hasNext();
		}

		@Override
		public VariantContext next() {
			return this.progress.watch(this.delegate.next());
		}

		@Override
		public void close()  {
			this.delegate.close();
			this.progress.finish();
			
		}

		//@Override public AbstractVCFCodec getCodec() { return this.delegate.getCodec(); }

		@Override public VCFHeader getHeader() { return this.delegate.getHeader(); }

		@Override
		public VariantContext peek() {
			return this.delegate.peek();
		}
    	
    	}
    
    public static void copyHeaderAndVariantsTo(final VCFIterator iter,final VariantContextWriter w) {
    	w.writeHeader(iter.getHeader());
    	copyVariantsTo(iter,w);
    	}

    
    public static long copyVariantsTo(final VCFIterator iter,final VariantContextWriter w) {
    	long n=0L;
    	while(iter.hasNext()) {
    		w.add(iter.next());
    		++n;
    	}
    	return n;
    }
    /**
     * copy InfoHeaderType.ALLELE from ctx from annotationVarian 
     */
    public static List<Object> copyAlternateAlleleAnnotations(
    		final List<Allele> destAlternateAlleles,
    		final VariantContext annotationVariant,
    		final String infoAttributeName
    		)
		{
    	final List<Object> attributes=new ArrayList<>(destAlternateAlleles.size());
    	final List<Allele> galts=annotationVariant.getAlternateAlleles();
		final List<Object> gatts = annotationVariant.getAttributeAsList(infoAttributeName);
		for(final Allele a:destAlternateAlleles)
			{
			Object found=null;
			//final int idx=gnomadCtx.getAlleleIndex(a);//non idx(REF)==0
			final int idx=galts.indexOf(a);
			
			if(idx>=0) {
				if(idx<gatts.size() && gatts.get(idx)!=null && !gatts.get(idx).equals(".")) {
					found= gatts.get(idx);
					}
				}
			attributes.add(found);
			}
		return attributes;
		}
	}
