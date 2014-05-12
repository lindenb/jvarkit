package com.github.lindenb.jvarkit.tools.vcfmerge;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.LineReaderUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterFactory;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import htsjdk.samtools.cmdline.Option;
import htsjdk.samtools.cmdline.StandardOptionDefinitions;
import htsjdk.samtools.cmdline.Usage;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SortingCollection;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryFactory;

public class VCFMerge extends AbstractCommandLineProgram
	{
	private static final Log LOG=Log.getInstance(VCFMerge.class);
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+" Merge VCFs";

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME,
    		doc="VCF files to process.",
    		minElements=1,
    		optional=false)
	public List<File> IN=new ArrayList<File>();

    @Option(shortName= StandardOptionDefinitions.REFERENCE_SHORT_NAME,
    		doc="Reference file.",
    		optional=false
    		)
	public File REF=null;

    
    @Option(shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME,
    		doc="VCF output. ",
    		optional=true)
	public File OUT=null;
    
	private SAMSequenceDictionary dictionary;
	
	private static class VCFHandler
		{
		VCFCodec vcfCodec = new VCFCodec();
		VCFHeader header=null;
		
		
		
		
		VariantContext parse(String line)
			{
			return vcfCodec.decode(line);
			}	
		}
	
	private List<VCFHandler> vcfHandlers=new ArrayList<VCFHandler>();
	
	private class VariantOfFile
		implements Comparable<VariantOfFile>
		{
		int fileIndex=-1;
		String line=null;
		private VariantContext var=null;
		
		boolean same(VariantOfFile var)
			{
			VariantContext vc1=parse();
			VariantContext vc2=var.parse();
			if(!vc1.getChr().equals(vc2.getChr())) return false;
			if(vc1.getStart()!=vc2.getStart()) return false;
			if(vc1.getEnd()!=vc2.getEnd()) return false;
			if(!vc1.getReference().equals(vc2.getReference())) return false;
			return true;
			}
		
		private int ref()
			{
			String chrom=parse().getChr();
			int refId=dictionary.getSequenceIndex(chrom);
			if(refId==-1) throw new RuntimeException("unknown chromosome "+ chrom+" in "+line);
			return refId;
			}

		
		public int compareTo(VariantOfFile var)
			{
			VariantContext vc1=parse();
			VariantContext vc2=var.parse();
			int i=ref()-var.ref();
			if(i!=0) return i;
			i=vc1.getStart()-vc2.getStart();
			if(i!=0) return i;
			i=vc1.getEnd()-vc2.getEnd();
			if(i!=0) return i;
			i=vc1.getReference().compareTo(vc2.getReference());
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

	private VariantContext buildContext(VCFCodec codec,VCFHeader header,List<VariantOfFile> row)
		{
		Double qual=null;
		Set<String> filters=new HashSet<String>();
		Set<Allele> alleles=new HashSet<Allele>();
		HashMap<String,Genotype> sample2genotype=new HashMap<String,Genotype>();
		VariantContextBuilder vcb=new VariantContextBuilder();
		Map<String,Object> atts=new HashMap<String,Object>();
		String id=null;
		vcb.chr(row.get(0).parse().getChr());
		vcb.start(row.get(0).parse().getStart());
		vcb.stop(row.get(0).parse().getEnd());
		
		for(VariantOfFile var:row)
			{
			VariantContext ctx=var.parse();
			
			if(qual==null || qual<ctx.getLog10PError()) qual=ctx.getLog10PError();
			
			filters.addAll(ctx.getFilters());
			alleles.addAll(ctx.getAlleles());
			if(id==null) id=ctx.getID();
			
			Map<String,Object> at=ctx.getAttributes();
			atts.putAll(at);
			
			for(String sample:ctx.getSampleNames())
				{
				Genotype g1=ctx.getGenotype(sample);
				if(g1==null) continue;
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
	protected int doWork()
		{
		PrintStream out=System.out;
		try {
	    	this.dictionary=new SAMSequenceDictionaryFactory().load(REF);
	
			Set<String> seenAttributes=new HashSet<String>();
			SortingCollection<VariantOfFile> array=SortingCollection.newInstance(
					VariantOfFile.class,
					new VariantCodec(),
					new VariantComparator(),
					super.MAX_RECORDS_IN_RAM
					);
			array.setDestructiveIteration(true);
			for(int fileIndex=0;fileIndex< this.IN.size();++fileIndex)
				{
				long nLine=0L;
				StringWriter sw=new StringWriter();
				File vcfFile= this.IN.get(fileIndex);
				LOG.info("reading from "+vcfFile+" "+( fileIndex+1)+"/"+ this.IN.size());
				VCFHandler handler=new VCFHandler();
				vcfHandlers.add(handler);
				
				BufferedReader in=IOUtils.openFileForBufferedReading(vcfFile);
				String line;
				while((line=in.readLine())!=null)
					{
					if(line.startsWith("#"))
						{
						sw.append(line).append('\n');
						continue;
						}
					if(handler.header==null)
						{
						ByteArrayInputStream bais=new ByteArrayInputStream(sw.toString().getBytes());
						LineIterator li=new LineIteratorImpl(LineReaderUtil.fromBufferedStream(bais));
						handler.header=(VCFHeader)handler.vcfCodec.readActualHeader(li);
						bais.close();
						}
					VariantOfFile vof=new VariantOfFile();
					vof.fileIndex=fileIndex;
					vof.line=line;
					VariantContext ctx=vof.parse();
					for(String att:ctx.getAttributes().keySet())
						{
						seenAttributes.add(att);
						}
					++nLine;
					if( nLine % 10000 ==0) LOG.info("reading var "+nLine +" of "+ vcfFile);
					array.add(vof);
					}
				
				in.close();
				}
			array.doneAdding();
			LOG.info("merging..."+vcfHandlers.size()+" vcfs");
			
			/* CREATE THE NEW VCH Header */
			VCFCodec mergeCodec=new VCFCodec();
			VCFHeader mergeHeader=null;
				{
				/* create the meta data */
				Set<VCFHeaderLine> metaData=new HashSet<VCFHeaderLine>();
				Set<String> sampleNames=new TreeSet<String>();
				for(VCFHandler handler: this.vcfHandlers)
					{
					metaData.addAll(handler.header.getMetaDataInSortedOrder());
					sampleNames.addAll(handler.header.getSampleNamesInOrder());
					for(VCFInfoHeaderLine vihl:handler.header.getInfoHeaderLines())
						{
						seenAttributes.remove(vihl.getKey());
						}
					}
				
				
				
				mergeHeader=new VCFHeader(
						metaData,
						sampleNames
						);
				
				//fix missing atts
				for(String att:seenAttributes)
					{
					mergeHeader.addMetaDataLine(new VCFInfoHeaderLine(att, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Undefined"));
					}
				
				}
			
			if(OUT!=null)
				{
				LOG.info("opening "+OUT);
				out=new PrintStream(IOUtils.openFileForWriting(OUT));
				}
				
			//create the context writer
			VariantContextWriter w= VariantContextWriterFactory.create(out,null,EnumSet.noneOf(Options.class));
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
					LOG.info("end of iteration");
					if(!row.isEmpty())
						{
						w.add( buildContext(mergeCodec,mergeHeader,row));
						}
					
					break;
					}
				
				if(!row.isEmpty() && !row.get(0).same(var))
					{
					w.add( buildContext(mergeCodec,mergeHeader,row));
					row.clear();
					}
				
				
				row.add(var);
					
				}
			w.close();
			LOG.info("done");
			}
		catch(Exception err)
			{
			LOG.error(err, ""+err.getMessage());
			return -1;
			}
		finally
			{
			out.flush();
			if(OUT!=null) out.close();
			}		
		return 0;
		}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new VCFMerge().instanceMainWithExit(args);

		}

	}
