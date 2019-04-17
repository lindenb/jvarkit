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
package com.github.lindenb.jvarkit.tools.vcfmerge;

import java.io.Closeable;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.LineReader;
import htsjdk.tribble.readers.SynchronousLineReader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.StringUtil;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;
/*
BEGIN_DOC

 
## Example


```bash
$  find ./ -name "*.vcf.gz" | xargs java -jar dist/vcfmerge.jar   > out.vcf
```

END_DOC
 */
@Program(name="vcfmerge",
	description="Merge VCF Files",
	deprecatedMsg="use GATK combineVariants.",
	keywords={"vcf","sort"}
	)
public class VCFMerge
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VCFMerge.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-s","--sorted"},description="files are known to be ROD sorted")
	private boolean filesAreSorted = false;

	@Parameter(names={"-m","--nomerge"},description="Do NOT merge VariantContext lines, but create multiple lines")
	private boolean doNotMergeRowLines = false;

	@Parameter(names={"-homref","--homref"},description="Use HomRef 0/0 for unknown variant")
	private boolean useHomRefForUnknown = false;
	
	@Parameter(names={"-region","--region"},description="Merge in that region: " + IntervalParser.OPT_DESC )
	private String regionStr = "";

	
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	
	
	
	private static final VCFInfoHeaderLine NO_MERGE_INFO_HEADER=
			new VCFInfoHeaderLine("VcfMergeOrigin",1,VCFHeaderLineType.String, "VCFmerge: origin of variant");
	

	/** user input files */
	private Set<String> userVcfFiles=new HashSet<String>();
	/** list<codec + header> */
	private List<VCFHandler> vcfHandlers=new ArrayList<VCFHandler>();
	
	

	public VCFMerge()
		{
		}
	
	
	
	private static class VCFHandler
		{
		final String origin;
		AbstractVCFCodec vcfCodec = VCFUtils.createDefaultVCFCodec();
		VCFHeader header=null;
		
		
		VCFHandler(final String origin) {
			this.origin=origin;
		}
		
		VariantContext parse(final String line)
			{
			return vcfCodec.decode(line);
			}
	
		
		
		}
	
	
	/** VariantContext associated to a file index */
	private class VariantOfFile
		implements Comparable<VariantOfFile>
		{
		/** will be used to retrieve VCFHeader */
		int fileIndex=-1;
		/** vcf line */
		String line=null;
		/** variantContext cache */
		private VariantContext var=null;
		
		boolean same(final VariantOfFile var)
			{
			return VCFMerge.this.compareChromPosRef.compare(this.parse(),var.parse())==0;
			}
		

		@Override
		public int compareTo(final VariantOfFile var)
			{
			final VariantContext vc1=parse();
			final VariantContext vc2=var.parse();
			final int i= VCFMerge.this.compareChromPosRef.compare(vc1, vc2);
			if(i!=0) return i;
			return fileIndex - var.fileIndex;
			}
		
		
		VariantContext parse()
			{
			if(this.var==null)
				{
				this.var=vcfHandlers.get(fileIndex).parse(this.line);
				}
			return var;
			}	
		}
	
	/**
	 * Variant serializer for sorting collection
	 *
	 */
	private  class VariantCodec
		extends AbstractDataCodec<VariantOfFile>
		{
		@Override
		public VariantOfFile decode(final DataInputStream dis) throws IOException
			{
			final VariantOfFile o=new VariantOfFile();
			try
				{
				o.fileIndex=dis.readInt();
				}
			catch(IOException err)
				{
				return null;
				}
			o.line=readString(dis);
			return o;

			}
		@Override
		public void encode(final DataOutputStream dos,final VariantOfFile s)
				throws IOException {
			dos.writeInt(s.fileIndex);
			writeString(dos,s.line);
			}
		@Override
		public VariantCodec clone() {
			return new VariantCodec();
			}
		}
	/** sorter for sorting collection */
	private class VariantComparator implements Comparator<VariantOfFile>
		{
		@Override
		public int compare(final VariantOfFile o1,final VariantOfFile o2)
			{
			return o1.compareTo(o2);
			}
		}
	
	private List<VariantContext> buildContextFromVariantOfFiles(
			VCFHeader header,
			List<VariantOfFile> row
			)
		{
		if(this.doNotMergeRowLines) {
			final List<VariantContext> L = new ArrayList<>(row.size());
			for(final VariantOfFile vof:row)
				{
				final VariantContext ctx = vof.parse();
				final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				final List<Genotype> genotypes = new ArrayList<>(ctx.getGenotypes());
				final Set<String> remainingSamples=new HashSet<String>(header.getSampleNamesInOrder());
				remainingSamples.removeAll(ctx.getSampleNames());
				for(String sampleName:remainingSamples)
					{
					genotypes.add(createMissingGenotype(sampleName,ctx.getReference()));
					}
				vcb.genotypes(genotypes);
				vcb.attribute(NO_MERGE_INFO_HEADER.getID(),
						VCFUtils.escapeInfoField(this.vcfHandlers.get(vof.fileIndex).origin)
						);
				L.add(vcb.make());
				}
			return L;
			}
		
		
		final List<VariantContext> row2=new ArrayList<VariantContext>(row.size());
		for(VariantOfFile vof:row) row2.add(vof.parse());
		return buildContextFromVariantContext(header,row2);
		}
	
	private Genotype createMissingGenotype(final String sampleName,final Allele ref)
		{
		if(this.useHomRefForUnknown)
			{
			return GenotypeBuilder.create(sampleName, Arrays.asList(
					ref,
					ref
					));
			}
		else
			{
			return GenotypeBuilder.createMissing(sampleName, 2);
			}		
		}
	
	private  Comparator<Genotype> genotypeComparator = (g1,g2) -> {
		 if(g2.hasGQ() && g1.hasGQ() )
		 	{
			return g2.getGQ() - g1.getGQ(); 
		 	}
		 return 0;
		};
		
	
	private List<VariantContext> buildContextFromVariantContext(
			final VCFHeader header,
			final List<VariantContext> row
			)
		{
		final VariantContextBuilder vcb=new VariantContextBuilder();
		final Map<String,Object> atts=new HashMap<String,Object>();
		final VariantContext ctx0 = row.get(0);
		
		vcb.chr(ctx0.getContig());
		vcb.start(ctx0.getStart());
		vcb.stop(ctx0.getEnd());
		
		//fill genotypes
		final HashMap<String,Genotype> sample2genotype=new HashMap<String,Genotype>();
		for(final VariantContext ctx:row)
			{
			for(final String sample:ctx.getSampleNames())
				{
				final Genotype g1=ctx.getGenotype(sample);
				if(g1==null || !g1.isCalled()) continue;
				final Genotype g2=sample2genotype.get(sample);
				if(g2==null || this.genotypeComparator.compare(g1, g2)<0)
					{
					sample2genotype.put(sample,g1);
					}
				}
			}
		// missing samples ?
		final Set<String> remainingSamples=new HashSet<String>(header.getSampleNamesInOrder());
		remainingSamples.removeAll(sample2genotype.keySet());
		for(final String sampleName : remainingSamples)
			{
			final Genotype missing = createMissingGenotype(sampleName,row.get(0).getReference());
			sample2genotype.put(sampleName,missing);
			}
		
		//collect alleles
		final List<Allele> alleleList =new ArrayList<>();
		alleleList.add(ctx0.getReference());
		for(final String sampleName:sample2genotype.keySet())
			{	
			final Genotype g1=sample2genotype.get(sampleName);
			if(!g1.isCalled() ) continue;
			for(final Allele ga: g1.getAlleles())
				{
				if(ga.isReference() || alleleList.contains(ga)) continue;
				alleleList.add(ga);
				}
			}
		
		
		
		vcb.attributes(atts);
		vcb.alleles(alleleList);
		vcb.genotypes(sample2genotype.values());
		return Collections.singletonList(vcb.make());
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		InputStream in=null;
		try
			{
			this.userVcfFiles.addAll(IOUtils.unrollFiles(args));
			
			if(this.userVcfFiles.isEmpty())
				{
				LOG.error("No input");
				return -1;
				}
			else if(this.userVcfFiles.size()==1)
				{
				in=IOUtils.openURIForReading(this.userVcfFiles.iterator().next());
				copyTo(in);
				in.close();
				in=null;
				}
			else
				{
				return workUsingPeekOrSorting();
				}
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			this.userVcfFiles.clear();
			}

		}
	
	private int workUsingPeekOrSorting() throws IOException
		{
		if(this.filesAreSorted)
			{
			return workUsingPeekIterator();
			}
		else
			{
			return workUsingSortingCollection();
			}
		}
	
	
	
	private void copyTo(final InputStream in) throws IOException
		{
		final VCFIterator iter= VCFUtils.createVCFIteratorFromInputStream(in);
		final VariantContextWriter out= this.openVariantContextWriter(outputFile);
		VCFUtils.copyHeaderAndVariantsTo(iter, out);
		CloserUtil.close(out);
		CloserUtil.close(iter);
		}
	
	
	/** container uri+vcfIterator */
	private class PeekVCF implements Closeable
		{
		final String uri;
		final VCFFileReader reader;
		final PeekableIterator<VariantContext> iter;
		final CloseableIterator<VariantContext> iter0;
		final VCFHeader header;
		final List<VariantContext> buffer = new ArrayList<>();
		
		PeekVCF(final String uri) throws IOException {
			this.uri = uri;
			if(StringUtil.isBlank(VCFMerge.this.regionStr))
				{
				this.reader = new VCFFileReader(new File(uri),false);
				this.header = this.reader.getFileHeader();
				this.iter0  = this.reader.iterator();
				}
			else
				{
				this.reader = new VCFFileReader(new File(uri),true);
				this.header = this.reader.getFileHeader();
				final IntervalParser intervalParser=new IntervalParser(this.header.getSequenceDictionary());
				intervalParser.setContigNameIsWholeContig(true);
				final Interval rgn = intervalParser.parse(VCFMerge.this.regionStr);
				this.iter0  = this.reader.query(rgn.getContig(), rgn.getStart(), rgn.getEnd());
				}
			this.iter = new PeekableIterator<>(this.iter0); 
			}
		
		private List<VariantContext> priv_peek()
			{
			if(!this.buffer.isEmpty()) return this.buffer;
			while(this.iter.hasNext())
				{
				final VariantContext ctx= this.iter.peek();
				if(this.buffer.isEmpty())
					{
					this.buffer.add(this.iter.next());
					}
				else
					{
					// compare with first item in buffer
					final int i = VCFMerge.this.compareChromPos.compare(
							ctx,
							this.buffer.get(0) 
							);
					if( i< 0) {
						throw new JvarkitException.UserError("Variant are not sorted! got: "+ctx+" after "+buffer.get(0));
						}
					else if(i > 0)
						{
						break;
						}
					else //i == 0 or growing
						{
						buffer.add(this.iter.next());
						}
					}
				}
			Collections.sort(this.buffer, VCFMerge.this.compareChromPosRef);
			return buffer;
			}
		List<VariantContext> peek()
			{
			final List<VariantContext> L = priv_peek();
			if(L.isEmpty() || L.size()==1) return L;
			return L.stream().
					filter(V->V==L.get(0) /* compare ptr */|| VCFMerge.this.compareChromPosRef.compare(L.get(0),V)==0).
					collect(Collectors.toList());
			}
		
		void reset(final VariantContext ctx0)
			{
			this.buffer.removeIf(V->
					 V.getContig().equals(ctx0.getContig()) &&
					 V.getStart() == ctx0.getStart()  &&
				     V.getReference().equals(ctx0.getReference()));
			}
		
		@Override
		public void close()
			{
			CloserUtil.close(this.iter);
			CloserUtil.close(this.iter0);
			CloserUtil.close(this.reader);
			}
		@Override
		public String toString() {
			return this.uri;
			}
		}
		
	
	private int workUsingPeekIterator()
		{
		VariantContextWriter out = null;
		final List<PeekVCF> input=new ArrayList<PeekVCF>();

		try {
			final Set<String> genotypeSampleNames=new TreeSet<String>();
			final Set<VCFHeaderLine> metaData=new HashSet<VCFHeaderLine>();
			
			//get all VCF, check same dict
			for(final String arg:this.userVcfFiles )
				{
				LOG.info("Opening "+arg);
				final PeekVCF p=new PeekVCF(arg);
				input.add(p);
				genotypeSampleNames.addAll(p.header.getSampleNamesInOrder());
				metaData.addAll(p.header.getMetaDataInInputOrder());
				if(this.global_dictionary==null)
					{
					this.global_dictionary= SequenceDictionaryUtils.extractRequired(p.header);
					}
				else if(!SequenceUtil.areSequenceDictionariesEqual(this.global_dictionary, SequenceDictionaryUtils.extractRequired(p.header)))
					{
					throw new JvarkitException.DictionariesAreNotTheSame(this.global_dictionary, p.header.getSequenceDictionary());
					}
				}
			
			if(this.global_dictionary==null)
				{
				throw new IllegalStateException("No Dict");
				}
			
			//super.addMetaData(metaData);
			if(this.doNotMergeRowLines)
				{
				metaData.add(NO_MERGE_INFO_HEADER);
				}
			final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(this.global_dictionary);
			out = super.openVariantContextWriter(this.outputFile);
			final VCFHeader headerOut=new VCFHeader(
					metaData,
					genotypeSampleNames
					);
			
			
			out.writeHeader(headerOut);
			final List<VariantContext> row=new ArrayList<VariantContext>(input.size());
			//find smallest ordered variant
			long nCountForGC=0L;
			for(;;)
				{
				row.clear();
				for(final PeekVCF peekVcf: input)
					{
					final List<VariantContext> peekVcfList = peekVcf.peek(); 
					if(peekVcf.peek().isEmpty()) continue;
					final VariantContext ctx= peekVcfList.get(0);
					if(row.isEmpty())
						{
						row.addAll(peekVcfList);
						continue;
						}
					final int cmp= this.compareChromPosRef.compare( ctx, row.get(0));
					if(cmp==0)
						{
						row.addAll(peekVcfList);
						}
					else if(cmp<0)
						{
						row.clear();
						row.addAll(peekVcfList);
						}
					}
				if(row.isEmpty()) break;
				
				for(final VariantContext merged: buildContextFromVariantContext(headerOut, row))
					{
					out.add(progress.watch(merged));
					}
				
				//consumme peeked variants
				for(final PeekVCF peekVcf: input)
					{
					peekVcf.reset(row.get(0));
					}
				if(nCountForGC++%100000==0) System.gc();
				}
			for(final PeekVCF peekVcf: input)
				{
				peekVcf.close();
				}
			input.clear();
			CloserUtil.close(out); out=null;
			progress.finish();
			
			LOG.info("Done peek sorting");
			return RETURN_OK;
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
		}
		finally
			{
			CloserUtil.close(out);
			for(final PeekVCF p: input)
				{
				p.close();
				}
			}
		}
	
	private SAMSequenceDictionary global_dictionary=null;
	
	private final Function<String,Integer> contig2tid=C->{
		final int tid = global_dictionary.getSequenceIndex(C);
		if(tid==-1) throw new JvarkitException.ContigNotFoundInDictionary(C, global_dictionary);
		return tid;
		};
	
	private final Comparator<String> compareContigs = (C1,C2)->{
		if(C1.equals(C2)) return 0;
		return contig2tid.apply(C1) - contig2tid.apply(C2);
		};
	
	private final Comparator<VariantContext> compareChromPos = (V1,V2)->{
		int i = compareContigs.compare(V1.getContig(),V2.getContig());
		if( i!=0 ) return i;
		return V1.getStart() - V2.getStart();
		};	
	private final Comparator<VariantContext> compareChromPosRef = (V1,V2)->{
		int i = compareChromPos.compare(V1,V2);
		if( i!=0 ) return i;
		return V1.getReference().compareTo(V2.getReference());
		};	
		
	
	private int workUsingSortingCollection() 
		{
		VariantContextWriter w=null;
		SortingCollection<VariantOfFile> array = null;
		InputStream in = null;
		CloseableIterator<VariantOfFile> iter=null;
			try {
			final List<String> IN=new ArrayList<String>(this.userVcfFiles);
			final Set<String> genotypeSampleNames=new TreeSet<String>();
			final Set<VCFHeaderLine> metaData=new HashSet<VCFHeaderLine>();
			array= SortingCollection.newInstance(
					VariantOfFile.class,
					new VariantCodec(),
					new VariantComparator(),
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			array.setDestructiveIteration(true);
			
			for(int fileIndex=0;fileIndex<  IN.size();++fileIndex)
				{
				final String vcfFile= IN.get(fileIndex);
				LOG.info("reading from "+vcfFile+" "+( fileIndex+1)+"/"+ IN.size());
				final VCFHandler handler=new VCFHandler(vcfFile);
				vcfHandlers.add(handler);
	
				in=IOUtils.openURIForReading(vcfFile);
				final LineReader lr=new SynchronousLineReader(in);
				final LineIterator lit=new LineIteratorImpl(lr);
				handler.header=(VCFHeader)handler.vcfCodec.readActualHeader(lit);
	
	
				final SAMSequenceDictionary dict1=handler.header.getSequenceDictionary();
				if(dict1==null) throw new RuntimeException("dictionary missing in "+vcfFile);
				if(dict1.isEmpty()) throw new RuntimeException("dictionary is Empty in "+vcfFile);
				genotypeSampleNames.addAll(handler.header.getSampleNamesInOrder());
				metaData.addAll(handler.header.getMetaDataInInputOrder());
				if(fileIndex==0)
					{
					this.global_dictionary=dict1;
					}
				else if(!SequenceUtil.areSequenceDictionariesEqual(global_dictionary, dict1))
					{
					throw new JvarkitException.DictionariesAreNotTheSame(global_dictionary, dict1);
					}
				final Predicate<VariantOfFile> accept;
				if(!StringUtil.isBlank(VCFMerge.this.regionStr)) {
					final IntervalParser intervalParser=new IntervalParser(dict1);
					intervalParser.setContigNameIsWholeContig(true);
					final Interval rgn = intervalParser.parse(VCFMerge.this.regionStr);
					accept = (VOL)->{
						final VariantContext ctx = VOL.parse();
						return rgn.intersects(new Interval(ctx.getContig(), ctx.getStart(),ctx.getEnd()));
						};
					}
				else
					{
					accept = (VOL) -> true;	
					}

				
				while(lit.hasNext())
					{					
					final VariantOfFile vof=new VariantOfFile();
					vof.fileIndex=fileIndex;
					vof.line=lit.next();
					if(!accept.test(vof)) continue;
					array.add(vof);
					}
	
				in.close();
				in=null;
				}
			array.doneAdding();
			LOG.info("merging..."+vcfHandlers.size()+" vcfs");
	
			/* CREATE THE NEW VCH Header */
			VCFHeader mergeHeader=null;
			//super.addMetaData(metaData);
			
			if(this.doNotMergeRowLines)
				{
				metaData.add(NO_MERGE_INFO_HEADER);
				}
			
			mergeHeader=new VCFHeader(
					metaData,
					genotypeSampleNames
					);
			
			
	
			//create the context writer
			w= super.openVariantContextWriter(outputFile);
			w.writeHeader(mergeHeader);
			iter= array.iterator();
			final List<VariantOfFile> row=new ArrayList<VariantOfFile>();
			for(;;)
				{
				VariantOfFile var=null;
				if(iter.hasNext())
					{
					var=iter.next();
					}
				else
					{
					LOG.info("end of iteration");
					if(!row.isEmpty())
						{
						for(final VariantContext merged:  buildContextFromVariantOfFiles(mergeHeader,row))
							{
							w.add(merged);
							}
						}
					break;
					}
	
				if(!row.isEmpty() && !row.get(0).same(var))
					{
					for(final VariantContext merged:  buildContextFromVariantOfFiles(mergeHeader,row))
						{
						w.add(merged);
						}
					row.clear();
					}
				row.add(var);
				}
			CloserUtil.close(w);w=null;
			array.cleanup();array=null;
			CloserUtil.close(iter);iter=null;
			LOG.info("done");
			return RETURN_OK;
			}
		catch(Exception err) {
				LOG.error(err);
				return -1;
			}
		finally
			{
			this.userVcfFiles.clear();
			CloserUtil.close(w);
			CloserUtil.close(in);
			CloserUtil.close(iter);
			if(array!=null) array.cleanup();
			}
		}

	public static void main(final String[] args)
		{
		new VCFMerge().instanceMainWithExit(args);
		}
	
	}
