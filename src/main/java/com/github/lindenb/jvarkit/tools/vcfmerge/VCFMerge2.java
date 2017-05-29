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
* 2015 adapted for knime
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.vcfmerge;

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
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SortingCollection;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
/*
BEGIN_DOC

 
## Example

```bash
$  find ./ -name "*.vcf.gz" | xargs java -jar dist/vcfmerge.jar   > out.vcf
```

END_DOC
 */
@Program(name="vcfmerge",description="Merge VCF Files",deprecatedMsg="use GATK combineVariants ")
public class VCFMerge2
	extends Launcher
	{

	private static final Logger LOG = Logger.build(VCFMerge2.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-s","--sorted"},description="files are known to be ROD sorted")
	private boolean filesAreSorted = false;

	@Parameter(names={"-m","--nomerge"},description="Do NOT merge VariantContext lines, but create multiple lines")
	private boolean doNotMergeRowLines = false;

	@Parameter(names={"-homref","--homref"},description="Use HomRef 0/0 for unknown variant")
	private boolean useHomRefForUnknown = false;

	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	
	
	
	private static final VCFInfoHeaderLine NO_MERGE_INFO_HEADER=
			new VCFInfoHeaderLine("VcfMergeOrigin",1,VCFHeaderLineType.String, "VCFmerge: origin of variant");
	

	/** user input files */
	private Set<String> userVcfFiles=new HashSet<String>();
	/** list<codec + header> */
	private List<VCFHandler> vcfHandlers=new ArrayList<VCFHandler>();
	
	
	/** moving to public for knime */
	public VCFMerge2()
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
		
		boolean same(VariantOfFile var)
			{
			return compare(global_dictionary,this.parse(),var.parse())==0;
			}
		

		
		public int compareTo(VariantOfFile var)
			{
			final VariantContext vc1=parse();
			final VariantContext vc2=var.parse();
			final int i= compare(global_dictionary, vc1, vc2);
			if(i!=0) return i;
			return fileIndex - var.fileIndex;
			}
		
		
		VariantContext parse()
			{
			if(var==null)
				{
				var=vcfHandlers.get(fileIndex).parse(this.line);
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
	
	private List<VariantContext> buildContextFromVariantContext(
			final VCFHeader header,
			final List<VariantContext> row
			)
		{
		
		
		Double qual=null;
		final Set<String> filters=new HashSet<String>();
		final HashMap<String,Genotype> sample2genotype=new HashMap<String,Genotype>();
		final VariantContextBuilder vcb=new VariantContextBuilder();
		final Map<String,Object> atts=new HashMap<String,Object>();
		String id=null;
		vcb.chr(row.get(0).getContig());
		vcb.start(row.get(0).getStart());
		vcb.stop(row.get(0).getEnd());
		
		final VCFInfoHeaderLine dpInfo= header.getInfoHeaderLine("DP");
		int total_dp=0;
		
		for(VariantContext ctx:row)
			{
			if(ctx.hasLog10PError())
				{
				if(qual==null || qual<ctx.getLog10PError()) qual=ctx.getLog10PError();
				}
			int dp=ctx.getAttributeAsInt("DP",-1);
			if(dp!=-1)
				{
				total_dp+=dp;
				}
			
			filters.addAll(ctx.getFilters());
			if(id==null && ctx.hasID()) id=ctx.getID();
			
			Map<String,Object> at=ctx.getAttributes();
			for(String attName:at.keySet())
				{
				if(row.size()!=1)
					{
					final VCFInfoHeaderLine headerInfo = header.getInfoHeaderLine(attName);
					if(headerInfo==null) continue;
					if(headerInfo.getCountType()==VCFHeaderLineCount.A) 
						{
						//TODO
						continue;
						}
					if(headerInfo.getCountType()==VCFHeaderLineCount.G) 
						{
						//TODO
						continue;
						}
					}
				atts.put(attName, at.get(attName));
				}
			
			
			for(final String sample:ctx.getSampleNames())
				{
				final Genotype g1=ctx.getGenotype(sample);
				if(g1==null || !g1.isCalled()) continue;
				final Genotype g2=sample2genotype.get(sample);
				if(g2==null || (g2.hasGQ() && g1.hasGQ() && g2.getGQ()<g1.getGQ()))
					{
					sample2genotype.put(sample,g1);
					}
				}
			}
		
		
		final Set<Allele> alleles=new HashSet<Allele>();
		alleles.add(row.get(0).getReference());
		for(final String sampleName:sample2genotype.keySet())
			{	
			final Genotype g1=sample2genotype.get(sampleName);
			if(!g1.isAvailable()) continue;
			alleles.addAll(g1.getAlleles());
			}
		// missing samples ?
		final Set<String> remainingSamples=new HashSet<String>(header.getSampleNamesInOrder());
		remainingSamples.removeAll(sample2genotype.keySet());
		for(String sampleName:remainingSamples)
			{
			final Genotype missing = createMissingGenotype(sampleName,row.get(0).getReference());
			sample2genotype.put(sampleName,missing);
			}
		
		if( atts.containsKey("DP") &&
			dpInfo!=null && dpInfo.getCount()==1 && dpInfo.getType()==VCFHeaderLineType.Integer)
			{
			atts.remove("DP");
			atts.put("DP", total_dp);
			}
		
		vcb.attributes(atts);
		filters.remove(VCFConstants.PASSES_FILTERS_v4);
		if(!filters.isEmpty()) vcb.filters(filters);
		vcb.alleles(alleles);
		if(qual!=null) vcb.log10PError(qual);
		vcb.genotypes(sample2genotype.values());
		if(id!=null) vcb.id(id);
		return Collections.singletonList(vcb.make());
		}
	
	
	@Override
	public int doWork(List<String> args) {
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
			return RETURN_OK;
			}
		catch(Exception err)
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
		if(filesAreSorted)
			{
			return workUsingPeekIterator();
			}
		else
			{
			return workUsingSortingCollection();
			}
		}
	
	
	
	
	private void copyTo(InputStream in) throws IOException
		{
		final VcfIterator iter= VCFUtils.createVcfIteratorFromInputStream(in);
		final VariantContextWriter out= this.openVariantContextWriter(outputFile);
		VCFUtils.copyHeaderAndVariantsTo(iter, out);
		CloserUtil.close(out);
		CloserUtil.close(iter);
		}
	/** container uri+vcfIterator */
	private static class PeekVCF
		{
		String uri;
		VcfIterator iter=null;
		}
	
	private int compare(
			final SAMSequenceDictionary dict,
			final VariantContext me,
			final VariantContext other)
		{
		int ref0=dict.getSequenceIndex(me.getContig());
		if(ref0<0) throw new IllegalArgumentException("unknown chromosome not in sequence dictionary: "+me.getContig());
		int ref1=dict.getSequenceIndex(other.getContig());
		if(ref1<0) throw new IllegalArgumentException("unknown chromosome not in sequence dictionary: "+other.getContig());
		
		int i=ref0-ref1;
		if(i!=0) return i;
		i=me.getStart()-other.getStart();
		if(i!=0) return i;
		/*i=me.getEnd()-other.getEnd();
		if(i!=0) return i; no we can have indel at the same place*/ 
		i=me.getReference().compareTo(other.getReference());
		if(i!=0) return i;

		return 0;
		}
	
	private int workUsingPeekIterator()
		{
		VariantContextWriter out = null;
		final List<PeekVCF> input=new ArrayList<PeekVCF>();

		try {
			SAMSequenceDictionary dict=null;
			final Set<String> genotypeSampleNames=new TreeSet<String>();
			final Set<VCFHeaderLine> metaData=new HashSet<VCFHeaderLine>();
			
			//get all VCF, check same dict
			for(String arg:this.userVcfFiles )
				{
				PeekVCF p=new PeekVCF();
				p.uri=arg;
				LOG.info("Opening "+p.uri);
				p.iter=VCFUtils.createVcfIterator(p.uri);
				input.add(p);
				final SAMSequenceDictionary dict1=p.iter.getHeader().getSequenceDictionary();
				if(dict1==null) throw new RuntimeException("dictionary missing in "+p.uri);
				if(dict1.isEmpty()) throw new RuntimeException("dictionary is Empty in "+p.uri);
				genotypeSampleNames.addAll(p.iter.getHeader().getSampleNamesInOrder());
				metaData.addAll(p.iter.getHeader().getMetaDataInInputOrder());
				if(dict==null)
					{
					dict=dict1;
					}
				else if(!SequenceUtil.areSequenceDictionariesEqual(dict, dict1))
					{
					throw new RuntimeException("Not the same Sequence dictionaries "+input.get(0).uri+" / "+p.uri);
					}
				}
			//super.addMetaData(metaData);
			if(this.doNotMergeRowLines)
				{
				metaData.add(NO_MERGE_INFO_HEADER);
				}
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
			out = super.openVariantContextWriter(outputFile);
			final VCFHeader headerOut=new VCFHeader(
					metaData,
					genotypeSampleNames);
			
			
			out.writeHeader(headerOut);
			boolean peeked[]=new boolean[input.size()];
			final List<VariantContext> row=new ArrayList<VariantContext>();
			//find smallest ordered variant
			for(;;)
				{
				row.clear();
				Arrays.fill(peeked, false);
				for(int i=0;i< input.size();++i)
					{
					if(!input.get(i).iter.hasNext()) continue;
					final VariantContext ctx=input.get(i).iter.peek();
					if(row.isEmpty())
						{
						Arrays.fill(peeked, false);
						peeked[i]=true;
						row.add(ctx);
						continue;
						}
					final int cmp=compare(dict, ctx, row.get(0));
					if(cmp==0)
						{
						peeked[i]=true;
						row.add(ctx);
						}
					else if(cmp<0)
						{
						Arrays.fill(peeked, false);
						peeked[i]=true;
						row.clear();
						row.add(ctx);
						}
					}
				if(row.isEmpty()) break;
				
				for(final VariantContext merged: buildContextFromVariantContext(headerOut, row))
					{
					out.add(progress.watch(merged));
					}
				
				//consumme peeked variants
				for(int i=0;i< input.size();++i)
					{
					if(peeked[i]) input.get(i).iter.next();
					}
				
				}
			CloserUtil.close(out); out=null;
			progress.finish();
			
			LOG.info("Done peek sorting");
			return RETURN_OK;
			}
		catch(Exception err) {
			LOG.error(err);
			return -1;
		}
		finally
			{
			CloserUtil.close(out);
			for(PeekVCF p: input)
				{
				CloserUtil.close(p.iter);
				}
			}
		}
	
	private SAMSequenceDictionary global_dictionary=null;
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
					this.writingSortingCollection.getTmpDirectories()
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
					global_dictionary=dict1;
					}
				else if(!SequenceUtil.areSequenceDictionariesEqual(global_dictionary, dict1))
					{
					throw new RuntimeException("Not the same Sequence dictionaries "+IN.get(0)+" / "+vcfFile);
					}
	
	
				while(lit.hasNext())
					{					
					VariantOfFile vof=new VariantOfFile();
					vof.fileIndex=fileIndex;
					vof.line=lit.next();
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

	

	

	/**
	 * main
	 */
	public static void main(String[] args)
		{
		new VCFMerge2().instanceMainWithExit(args);
		}
	
	}
