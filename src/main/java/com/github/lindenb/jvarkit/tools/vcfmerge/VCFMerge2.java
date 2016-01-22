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

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.LineReader;
import htsjdk.tribble.readers.LineReaderUtil;
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
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SortingCollection;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.knime.KnimeApplication;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.VcfIteratorImpl;

/** merge VCF files */
public class VCFMerge2
	extends com.github.lindenb.jvarkit.util.AbstractCommandLineProgram
	implements KnimeApplication
	{
	/** output file : stdout if null */
	private File outputFile=null;
	/** number of variants filterer */
	private int countMergedVariants=0;
	/** input file are sorted */
	private boolean filesAreSorted=false;
	/** user input files */
	private Set<String> userVcfFiles=new HashSet<String>();
	/** list<codec + header> */
	private List<VCFHandler> vcfHandlers=new ArrayList<VCFHandler>();
	
	
	/** moving to public for knime */
	public VCFMerge2()
		{
		}
	
	
	@Override
	public void setOutputFile(File out) {
		this.outputFile=out;
		}
	
	public File getOutputFile() {
		return outputFile;
		}

	public void setFilesAreSorted(boolean filesAreSorted) {
		this.filesAreSorted = filesAreSorted;
		}
	
	@Override
	public String getProgramDescription()
		{
		return " Merge VCF Files.";
		}

    @Override
    protected String getOnlineDocUrl() {
    	return DEFAULT_WIKI_PREFIX+"VCFMerge";
    	}
    
    public int getCountMergedVariants() {
		return countMergedVariants;
		}
  
	
	private static class VCFHandler
		{
		AbstractVCFCodec vcfCodec = VCFUtils.createDefaultVCFCodec();
		VCFHeader header=null;
		

		VariantContext parse(String line)
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
			VariantContext vc1=parse();
			VariantContext vc2=var.parse();
			int i= compare(global_dictionary, vc1, vc2);
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
	
	private  class VariantCodec
		extends AbstractDataCodec<VariantOfFile>
		{
		@Override
		public VariantOfFile decode(DataInputStream dis) throws IOException
			{
			VariantOfFile o=new VariantOfFile();
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
		public void encode(DataOutputStream dos, VariantOfFile s)
				throws IOException {
			dos.writeInt(s.fileIndex);
			writeString(dos,s.line);
			}
		@Override
		public VariantCodec clone() {
			return new VariantCodec();
			}
		}
	
	private class VariantComparator implements Comparator<VariantOfFile>
		{
		@Override
		public int compare(VariantOfFile o1, VariantOfFile o2)
			{
			return o1.compareTo(o2);
			}
		}
	
	private VariantContext buildContextFromVariantOfFiles(
			VCFHeader header,
			List<VariantOfFile> row
			)
		{
		List<VariantContext> row2=new ArrayList<VariantContext>(row.size());
		for(VariantOfFile vof:row) row2.add(vof.parse());
		return buildContextFromVariantContext(header,row2);
		}
	
	private VariantContext buildContextFromVariantContext(
			VCFHeader header,
			List<VariantContext> row
			)
		{
		Double qual=null;
		Set<String> filters=new HashSet<String>();
		Set<Allele> alleles=new HashSet<Allele>();
		HashMap<String,Genotype> sample2genotype=new HashMap<String,Genotype>();
		VariantContextBuilder vcb=new VariantContextBuilder();
		Map<String,Object> atts=new HashMap<String,Object>();
		String id=null;
		vcb.chr(row.get(0).getContig());
		vcb.start(row.get(0).getStart());
		vcb.stop(row.get(0).getEnd());
		
		VCFInfoHeaderLine dpInfo= header.getInfoHeaderLine("DP");
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
			alleles.addAll(ctx.getAlleles());
			if(id==null && ctx.hasID()) id=ctx.getID();
			
			Map<String,Object> at=ctx.getAttributes();
			atts.putAll(at);
			
			for(String sample:ctx.getSampleNames())
				{
				Genotype g1=ctx.getGenotype(sample);
				if(g1==null || !g1.isCalled()) continue;
				Genotype g2=sample2genotype.get(sample);
				if(g2==null || (g2.hasGQ() && g1.hasGQ() && g2.getGQ()<g1.getGQ()))
					{
					sample2genotype.put(sample,g1);
					}
				}
			}
		
		// missing samples ?
		Set<String> remainingSamples=new HashSet<String>(header.getSampleNamesInOrder());
		remainingSamples.removeAll(sample2genotype.keySet());
		for(String sampleName:remainingSamples)
			{
			sample2genotype.put(sampleName, GenotypeBuilder.createMissing(sampleName, 2));
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
		return vcb.make();
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -s files are known to be sorted");		
		out.println(" -f (list) read vcf path from file");	
		out.println(" -o (output file) default: stdout");
		super.printOptions(out);
		}

	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		
		while((c=opt.getopt(args,getGetOptDefault()+ "sf:o:"))!=-1)
			{
			switch(c)
				{
				case 'o': setOutputFile(new File(opt.getOptArg())); break;
				case 's': setFilesAreSorted(true);break;
				case 'f':
					{
					try {
						BufferedReader in=IOUtils.openURIForBufferedReading(opt.getOptArg());
						String line;
						while((line=in.readLine())!=null)
							{
							if(line.startsWith("#") || line.trim().isEmpty()) continue;
							this.userVcfFiles.add(line.trim());
							}
						in.close();
						} 
					catch (Exception e)
						{
						error(e);
						return -1;
						}
					break;
					}
				default: 
					{
					switch(handleOtherOptions(c, opt, args))
						{
						case EXIT_FAILURE:return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		try
			{
			if(opt.getOptInd()==args.length)
				{
				if(!userVcfFiles.isEmpty())
					{
					this.workUsingPeekOrSorting();
					}
				else
					{
					copyTo(System.in);
					}
				}
			else if(opt.getOptInd()+1==args.length)
				{
				if(!this.userVcfFiles.isEmpty())
					{
					error("option -f used but filenames provided too.");
					error(getMessageBundle("illegal.number.of.arguments"));
					return -1;
					}
				String filename=args[opt.getOptInd()];
				InputStream in=IOUtils.openURIForReading(filename);
				copyTo(in);
				CloserUtil.close(in);
				}
			else
				{
				for(int i=opt.getOptInd();i< args.length;++i)
					{
					this.userVcfFiles.add(args[i]);
					}
				
				this.workUsingPeekOrSorting();
				}
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		

		}
	
	private void workUsingPeekOrSorting() throws IOException
		{
		if(filesAreSorted)
			{
			workUsingPeekIterator();
			}
		else
			{
			workUsingSortingCollection();
			}
		}
	
	private VariantContextWriter createVariantContextWriter()
		throws IOException
		{
		return this.outputFile==null?
				VCFUtils.createVariantContextWriterToStdout():
				VCFUtils.createVariantContextWriter(this.outputFile)
				;
		}
	
	
	private void copyTo(InputStream in) throws IOException
		{
		this.countMergedVariants=0;
		@SuppressWarnings("resource")
		VcfIterator iter=new VcfIteratorImpl(in);
		VCFHeader h=iter.getHeader();
		VariantContextWriter out= this.createVariantContextWriter();
		out.writeHeader(h);
		while(iter.hasNext())
			{
			this.countMergedVariants++;
			out.add(iter.next());
			}
		CloserUtil.close(out);
		}
	
	private static class PeekVCF
		{
		String uri;
		VcfIterator iter=null;
		}
	private int compare(
			SAMSequenceDictionary dict,
			VariantContext me,
			VariantContext other)
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
	
	private void workUsingPeekIterator()
		throws IOException
		{
		this.countMergedVariants=0;
		SAMSequenceDictionary dict=null;
		Set<String> genotypeSampleNames=new HashSet<String>();
		Set<VCFHeaderLine> metaData=new HashSet<VCFHeaderLine>();
		List<PeekVCF> input=new ArrayList<PeekVCF>();
		
		//get all VCF, check same dict
		for(String arg:this.userVcfFiles )
			{
			PeekVCF p=new PeekVCF();
			p.uri=arg;
			info("Opening "+p.uri);
			p.iter=VCFUtils.createVcfIterator(p.uri);
			input.add(p);
			SAMSequenceDictionary dict1=p.iter.getHeader().getSequenceDictionary();
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
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
		VariantContextWriter out= this.createVariantContextWriter();
		VCFHeader headerOut=new VCFHeader(
				metaData,
				genotypeSampleNames);
		
		headerOut.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		headerOut.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		headerOut.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
		headerOut.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));

		
		out.writeHeader(headerOut);
		boolean peeked[]=new boolean[input.size()];
		List<VariantContext> row=new ArrayList<VariantContext>();
		//find smallest ordered variant
		for(;;)
			{
			row.clear();
			Arrays.fill(peeked, false);
			for(int i=0;i< input.size();++i)
				{
				if(!input.get(i).iter.hasNext()) continue;
				VariantContext ctx=input.get(i).iter.peek();
				if(row.isEmpty())
					{
					Arrays.fill(peeked, false);
					peeked[i]=true;
					row.add(ctx);
					continue;
					}
				int cmp=compare(dict, ctx, row.get(0));
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
			
			VariantContext merged= progress.watch( buildContextFromVariantContext(headerOut, row));
			out.add(merged);
			this.countMergedVariants++;

			//consumme peeked variants
			for(int i=0;i< input.size();++i)
				{
				if(peeked[i]) input.get(i).iter.next();
				}
			
			}
		CloserUtil.close(out);
		progress.finish();
		for(PeekVCF p: input)
			{
			CloserUtil.close(p.iter);
			}
		info("Done peek sorting");
		}
	
	private SAMSequenceDictionary global_dictionary=null;
	private int workUsingSortingCollection()
		throws IOException
		{
		List<String> IN=new ArrayList<String>(this.userVcfFiles);
		Set<String> genotypeSampleNames=new HashSet<String>();
		Set<VCFHeaderLine> metaData=new HashSet<VCFHeaderLine>();
		this.countMergedVariants=0;
		SortingCollection<VariantOfFile> array= SortingCollection.newInstance(
				VariantOfFile.class,
				new VariantCodec(),
				new VariantComparator(),
				10000,
				this.getTmpDirectories()
				);
		array.setDestructiveIteration(true);
		for(int fileIndex=0;fileIndex<  IN.size();++fileIndex)
			{
			String vcfFile= IN.get(fileIndex);
			info("reading from "+vcfFile+" "+( fileIndex+1)+"/"+ IN.size());
			VCFHandler handler=new VCFHandler();
			vcfHandlers.add(handler);

			InputStream in=IOUtils.openURIForReading(vcfFile);
			LineReader lr=LineReaderUtil.fromBufferedStream(in);
			LineIterator lit=new LineIteratorImpl(lr);
			handler.header=(VCFHeader)handler.vcfCodec.readActualHeader(lit);


			SAMSequenceDictionary dict1=handler.header.getSequenceDictionary();
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
		}
		array.doneAdding();
		info("merging..."+vcfHandlers.size()+" vcfs");

		/* CREATE THE NEW VCH Header */
		VCFHeader mergeHeader=null;

		mergeHeader=new VCFHeader(
				metaData,
				genotypeSampleNames
				);

		mergeHeader.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		mergeHeader.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		mergeHeader.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
		mergeHeader.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));


		//create the context writer
		VariantContextWriter w= this.createVariantContextWriter();
		w.writeHeader(mergeHeader);
		CloseableIterator<VariantOfFile> iter= array.iterator();
		List<VariantOfFile> row=new ArrayList<VariantOfFile>();
		for(;;)
			{
			VariantOfFile var=null;
			if(iter.hasNext())
				{
				var=iter.next();
				}
			else
				{
				info("end of iteration");
				if(!row.isEmpty())
					{
					this.countMergedVariants++;
					w.add( buildContextFromVariantOfFiles(mergeHeader,row));
					}
	
				break;
				}

			if(!row.isEmpty() && !row.get(0).same(var))
				{
				w.add( buildContextFromVariantOfFiles(mergeHeader,row));
				row.clear();
				}
			row.add(var);
			this.countMergedVariants++;
			}
		w.close();
		info("done");

		return 0;
		}

	
	@Override
	public int initializeKnime() {
		return 0;
	}

	@Override
	public int executeKnime(List<String> args)
		{
		this.userVcfFiles=new HashSet<String>(args);
		InputStream in=null;
		try
			{
			if(this.userVcfFiles.isEmpty())
				{
				error("No input");
				return -1;
				}
			else if(this.userVcfFiles.size()==1)
				{
				in=IOUtils.openURIForReading(args.get(0));
				copyTo(in);
				in.close();
				in=null;
				}
			else
				{
				workUsingPeekOrSorting();
				}
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			}
		return 0;
		}

	@Override
	public void disposeKnime() {
		
	}

	@Override
	public void checkKnimeCancelled() {
		
	}

	/**
	 * main
	 */
	public static void main(String[] args)
		{
		new VCFMerge2().instanceMainWithExit(args);
		}
	
	}
