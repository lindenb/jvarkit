package com.github.lindenb.jvarkit.tools.vcfmerge;

import java.io.DataInputStream;
import java.io.DataOutputStream;
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
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.util.picard.PicardException;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SortingCollection;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SortingCollectionFactory;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class VCFMerge2
	extends com.github.lindenb.jvarkit.util.AbstractCommandLineProgram
	{
	private SortingCollectionFactory<VariantOfFile> sortingCollectionFactory=new SortingCollectionFactory<VariantOfFile>();
	private List<VCFHandler> vcfHandlers=new ArrayList<VCFHandler>();

	
	private VCFMerge2()
		{
		
		}
	
	
	@Override
	public String getProgramDescription()
		{
		return " Merge VCF Files.";
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
	
	
	private class VariantOfFile
		implements Comparable<VariantOfFile>
		{
		int fileIndex=-1;
		String line=null;
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
			try
				{
				VariantOfFile o=new VariantOfFile();
				o.fileIndex=dis.readInt();
				o.line=dis.readUTF();
				return o;
				}
			catch(IOException err)
				{
				return null;
				}
			}
		@Override
		public void encode(DataOutputStream dos, VariantOfFile s)
				throws IOException {
			dos.writeInt(s.fileIndex);
			dos.writeUTF(s.line);
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
		if(row.size()==1) return row.get(0);
		Double qual=null;
		Set<String> filters=new HashSet<String>();
		Set<Allele> alleles=new HashSet<Allele>();
		HashMap<String,Genotype> sample2genotype=new HashMap<String,Genotype>();
		VariantContextBuilder vcb=new VariantContextBuilder();
		Map<String,Object> atts=new HashMap<String,Object>();
		String id=null;
		vcb.chr(row.get(0).getChr());
		vcb.start(row.get(0).getStart());
		vcb.stop(row.get(0).getEnd());
		
		for(VariantContext ctx:row)
			{
			if(ctx.hasLog10PError())
				{
				if(qual==null || qual<ctx.getLog10PError()) qual=ctx.getLog10PError();
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
		super.printOptions(out);
		}

	@Override
	public int doWork(String[] args)
		{
		boolean filesAreSorted=false;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "s"))!=-1)
			{
			switch(c)
				{
				case 's': filesAreSorted=true;break;
				default: 
					{
					switch(handleOtherOptions(c, opt, null))
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
				copyTo(System.in);
				return 0;
				}
			else if(opt.getOptInd()+1==args.length)
				{
				String filename=args[opt.getOptInd()];
				InputStream in=IOUtils.openURIForReading(filename);
				copyTo(in);
				CloserUtil.close(in);
				}
			else if(filesAreSorted)
				{
				workUsingPeekIterator(opt.getOptInd(),args);
				}
			else
				{
				workUsingSortingCollection(opt.getOptInd(),args);
				}
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		

		}
	
	
	
	private void copyTo(InputStream in) throws IOException
		{
		@SuppressWarnings("resource")
		VcfIterator iter=new VcfIterator(in);
		VCFHeader h=iter.getHeader();
		VariantContextWriter out=VCFUtils.createVariantContextWriter(null);
		out.writeHeader(h);
		while(iter.hasNext()) out.add(iter.next());
		CloserUtil.close(out);
		}
	
	private static class PeekVCF
		{
		String uri;
		InputStream is=null;
		VcfIterator iter=null;
		}
	private int compare(
			SAMSequenceDictionary dict,
			VariantContext me,
			VariantContext other)
		{
		int ref0=dict.getSequenceIndex(me.getChr());
		if(ref0<0) throw new IllegalArgumentException("unknown chromosome not in sequence dictionary: "+me.getChr());
		int ref1=dict.getSequenceIndex(other.getChr());
		if(ref1<0) throw new IllegalArgumentException("unknown chromosome not in sequence dictionary: "+other.getChr());
		
		int i=ref0-ref1;
		if(i!=0) return i;
		i=me.getStart()-other.getStart();
		if(i!=0) return i;
		i=me.getEnd()-other.getEnd();
		if(i!=0) return i;
		i=me.getReference().compareTo(other.getReference());
		if(i!=0) return i;

		return 0;
		}
	
	protected void workUsingPeekIterator(int optind,String args[])
		throws IOException
		{
		SAMSequenceDictionary dict=null;
		Set<String> genotypeSampleNames=new HashSet<String>();
		Set<VCFHeaderLine> metaData=new HashSet<VCFHeaderLine>();
		List<PeekVCF> input=new ArrayList<PeekVCF>();
		for(int i=optind;i< args.length;++i)
			{
			PeekVCF p=new PeekVCF();
			p.uri=args[i];
			info("Opening "+p.uri);
			p.is=IOUtils.openURIForReading(p.uri);
			p.iter=new VcfIterator(p.is);
			input.add(p);
			SAMSequenceDictionary dict1=p.iter.getHeader().getSequenceDictionary();
			if(dict1==null) throw new PicardException("dictionary missing in "+p.uri);
			if(dict1.isEmpty()) throw new PicardException("dictionary is Empty in "+p.uri);
			genotypeSampleNames.addAll(p.iter.getHeader().getSampleNamesInOrder());
			metaData.addAll(p.iter.getHeader().getMetaDataInInputOrder());
			if(i==optind)
				{
				dict=dict1;
				}
			else if(!SequenceUtil.areSequenceDictionariesEqual(dict, dict1))
				{
				throw new PicardException("Not the same Sequence dictionaries "+input.get(0).uri+" / "+p.uri);
				}
			}
		VariantContextWriter out=VCFUtils.createVariantContextWriter(null);
		
		VCFHeader headerOut=new VCFHeader(
				metaData,
				genotypeSampleNames);
		out.writeHeader(headerOut);
		boolean peeked[]=new boolean[input.size()];
		List<VariantContext> row=new ArrayList<VariantContext>();
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
			
			VariantContext merged=buildContextFromVariantContext(headerOut, row);
			out.add(merged);
			
			for(int i=0;i< input.size();++i)
				{
				if(peeked[i]) input.get(i).iter.next();
				}
			
			}
		CloserUtil.close(out);
		
		for(PeekVCF p: input)
			{
			CloserUtil.close(p.iter);
			CloserUtil.close(p.is);
			}
		}
	
	private SAMSequenceDictionary global_dictionary=null;
	protected int workUsingSortingCollection(int optind,String args[])
		throws IOException
		{
		List<String> IN=new ArrayList<String>();
		for(int i=optind;i< args.length;++i)
			{
			IN.add(args[i]);
			}
		this.sortingCollectionFactory.setComponentType(VariantOfFile.class);
		this.sortingCollectionFactory.setCodec(new VariantCodec());
		this.sortingCollectionFactory.setComparator(new VariantComparator());
		this.sortingCollectionFactory.setTmpDirs(this.getTmpDirectories());
		Set<String> genotypeSampleNames=new HashSet<String>();
		Set<VCFHeaderLine> metaData=new HashSet<VCFHeaderLine>();
	
		
		SortingCollection<VariantOfFile> array=this.sortingCollectionFactory.make();
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
				if(dict1==null) throw new PicardException("dictionary missing in "+vcfFile);
				if(dict1.isEmpty()) throw new PicardException("dictionary is Empty in "+vcfFile);
				genotypeSampleNames.addAll(handler.header.getSampleNamesInOrder());
				metaData.addAll(handler.header.getMetaDataInInputOrder());
				if(fileIndex==0)
					{
					global_dictionary=dict1;
					}
				else if(!SequenceUtil.areSequenceDictionariesEqual(global_dictionary, dict1))
					{
					throw new PicardException("Not the same Sequence dictionaries "+IN.get(0)+" / "+vcfFile);
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
				
				
				
			//create the context writer
			VariantContextWriter w= VCFUtils.createVariantContextWriter(null);
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
					
				}
			w.close();
			info("done");
					
		return 0;
		}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new VCFMerge2().instanceMainWithExit(args);

		}

	}
