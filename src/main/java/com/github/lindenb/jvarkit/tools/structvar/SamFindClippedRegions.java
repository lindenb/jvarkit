package com.github.lindenb.jvarkit.tools.structvar;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.broad.tribble.readers.LineIterator;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFConstants;
import org.broadinstitute.variant.vcf.VCFFormatHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFHeaderVersion;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.variant.vcf.VCFSimpleHeaderLine;


import net.sf.picard.PicardException;
import net.sf.picard.sam.MergingSamRecordIterator;
import net.sf.picard.sam.SamFileHeaderMerger;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.util.CloserUtil;
import net.sf.samtools.util.SequenceUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.Predicate;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SamSequenceRecordTreeMap;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

public class SamFindClippedRegions extends AbstractCommandLineProgram
	{
	private int min_clip_length=20;
	//private boolean ignore_poly_x=false;
	//private int rounding_pos=5;
	//private float min_fraction_of_clipped_records=0.3f;
	private int min_depth=15;
	

	
	private static class FilteringIterator<T>
		implements Iterator<T>
		{
		private Iterator<T> delegate;
		private Predicate<T> filter;
		private T _next=null;
		private boolean _hasNextCalled=false;
		FilteringIterator(Iterator<T> delegate,Predicate<T> filter)
			{
			this.delegate=delegate;
			this.filter=filter;
			}
		@Override
		public boolean hasNext()
			{
			if(_hasNextCalled) return _next!=null;
			_hasNextCalled=true;
			while(delegate.hasNext())
				{
				T o=delegate.next();
				if(filter.apply(o))
					{
					_next=o;
					break;
					}
				}
			return _next!=null;
			}
		@Override
		public T next()
			{
			if(!hasNext()) throw new IllegalStateException();
			_hasNextCalled=false;
			T o=_next;
			_next=null;
			return o;
			}
		@Override
		public void remove()
			{
			throw new UnsupportedOperationException();
			}
		}

	
	
	
	private class Input
		implements Closeable
		{
		//int index;
		File bamFile;
		SAMFileReader samFileReaderScan;
		SAMFileHeader header;
		String sampleName=null;
		Input(File bamFile)
			{
			info("Reading from "+bamFile);
			this.bamFile=bamFile;
			this.samFileReaderScan=new SAMFileReader(this.bamFile);
			this.samFileReaderScan.setValidationStringency(ValidationStringency.SILENT);
			this.header=this.samFileReaderScan.getFileHeader();
			for(SAMReadGroupRecord g:this.header.getReadGroups())
				{
				if(g.getSample()==null || g.getSample().isEmpty())
					{
					warning("Read group "+g.getId()+" has no sample in "+bamFile);
					continue;
					}
				if(sampleName!=null && !sampleName.equals(g.getSample()))
					{	
					throw new PicardException("Sorry one bam can only contain one sample and I got "+sampleName+" and "+g.getSample()+" in "+bamFile);
					}
				sampleName=g.getSample();
				}
			if(sampleName==null)
				{
				throw new PicardException("Bam "+this.bamFile+" doesn't contain a sample defined in header/read-group");
				}
			info("In "+bamFile+" Sample is "+sampleName);
			}
		
		@Override
		public void close() throws IOException {
			CloserUtil.close(samFileReaderScan);
			}
		}
	
	@Override
	public String getProgramDescription() {
		return "";
		}
	

	/*private static boolean closeTo(int pos1,int pos2, int max)
		{
		return Math.abs(pos2-pos1)<=max;
		}*/
	
	/*
	private static boolean same(char c1,char c2)
		{
		if(c1=='N' || c2=='N') return false;
		return Character.toUpperCase(c1)==Character.toUpperCase(c2);
		}*/
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -c (int) min size of clipped read. default:"+min_clip_length);
		out.println(" -x ignore poly-X");
		out.println(" -B (file) bed file (optional)");
		out.println(" -R (ref) "+getMessageBundle("reference.faidx")+". Optional");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		//final int CLOSE_DISTANCE=5;
		int readLength=150;
		File bedFile=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"c:B:"))!=-1)
			{
			switch(c)
				{
				case 'B': bedFile=new File(opt.getOptArg());break;
				//case 'x': ignore_poly_x=true;break;
				case 'c':min_clip_length=Integer.parseInt(opt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
	
		
		if(opt.getOptInd()==args.length)
			{
			error(getMessageBundle("illegal.number.of.arguments"));
			return -1;
			}
		
		List<Input> inputs=new ArrayList<Input>();
		
		VariantContextWriter w=null;
		//SAMFileWriter w=null;
		try
			{
			SAMSequenceDictionary dict=null;

			/* create input, collect sample names */
			Map<String,Input> sample2input=new HashMap<String,Input>();
			for(int optind=opt.getOptInd();
					optind< args.length;//TODO
					++optind)
				{
				Input input=new Input(new File(args[optind]));
				//input.index=inputs.size();
				inputs.add(input);
				
				if(sample2input.containsKey(input.sampleName))
					{
					error("Duplicate sample "+input.sampleName+" in "+input.bamFile+" and "+sample2input.get(input.sampleName).bamFile);
					return -1;
					}
				sample2input.put(input.sampleName, input);
				if(dict==null)
					{
					dict=input.header.getSequenceDictionary();
					}
				else if(!SequenceUtil.areSequenceDictionariesEqual(dict, input.header.getSequenceDictionary()))
					{
					error("Found more than one dictint sequence dictionary");
					return -1;
					}
					
				}
			info("Sample N= "+sample2input.size());
			
			
			/* create merged iterator */
			List<SAMFileHeader> headers=new ArrayList<SAMFileHeader>(sample2input.size());
			for(Input input:inputs) headers.add(input.header);
			SamFileHeaderMerger headerMerger=new SamFileHeaderMerger(SortOrder.coordinate,headers, true);
			
			List<SAMFileReader> readers=new ArrayList<SAMFileReader>(sample2input.size());
			for(Input input:inputs) readers.add(input.samFileReaderScan);
			
			MergingSamRecordIterator merginIter=new MergingSamRecordIterator(headerMerger, readers, true);
			
			Allele reference_allele=Allele.create("N",true);			
			Allele alternate_alleles[]=new Allele[]{ 
				Allele.create("<CLIP5>",false),
				Allele.create("<CLIP3>",false)
				};

			Set<VCFHeaderLine> vcfHeaderLines=new HashSet<VCFHeaderLine>();
			for(Allele alt:alternate_alleles)
				{
				vcfHeaderLines.add(new VCFSimpleHeaderLine(
						"<ID="+alt.getDisplayString()+",Description=\"StructVar\">",
						VCFHeaderVersion.VCF4_1,
						VCFConstants.ALT_HEADER_START.substring(2),
						Arrays.asList("ID", "Description")
						));
				}
			
			
			vcfHeaderLines.add(new VCFInfoHeaderLine("COUNT_SAMPLES", 1, VCFHeaderLineType.Integer, "Number of samples with  depth>="+this.min_depth));
			vcfHeaderLines.add(new VCFInfoHeaderLine(VCFConstants.DEPTH_KEY, 1, VCFHeaderLineType.Integer, "Approximate read depth."));
			vcfHeaderLines.add(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1, VCFHeaderLineType.String, "Genotype"));
			vcfHeaderLines.add(new VCFFormatHeaderLine(VCFConstants.DEPTH_KEY, 1, VCFHeaderLineType.Integer, "Approximate read depth"));
			vcfHeaderLines.add(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			vcfHeaderLines.add(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			for(int side=0;side<2;++side)
				{
				vcfHeaderLines.add(new VCFFormatHeaderLine("CN"+(side==0?5:3),1,VCFHeaderLineType.Integer,"count clipped in "+(side==0?5:3)+"'"));
				}
			if(dict!=null)
				{
				vcfHeaderLines.addAll(VCFUtils.samSequenceDictToVCFContigHeaderLine(dict));
				}
			
			VCFHeader vcfHeader=new VCFHeader(vcfHeaderLines,sample2input.keySet());
			
			w=VCFUtils.createVariantContextWriterToStdout();
			w.writeHeader(vcfHeader);
			

			
			SamSequenceRecordTreeMap<Boolean> intervals=null;
			//w=swf.make(header, System.out);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
			
			if(bedFile!=null)
				{
				intervals=new SamSequenceRecordTreeMap<Boolean>(dict);
				info("Reading "+bedFile);
				LineIterator r=IOUtils.openFileForLineIterator(bedFile);
				while(r.hasNext())
					{
					String line=r.next();
					if(line.startsWith("#") || line.isEmpty()) continue;
					String tokens[]=line.split("[\t]");
					if(tokens.length<3)
						{
						warning("Bad bed line in "+bedFile+" "+line);
						continue;
						}
					if(dict!=null && dict.getSequence(tokens[0])==null)
						{
						warning("undefined chromosome  in "+bedFile+" "+line);
						continue;
						}
					intervals.put(tokens[0],Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]), true);
					}
				CloserUtil.close(r);
				}
			
			
			
			LinkedList<SAMRecord> buffer=new LinkedList<SAMRecord>();
			
			final Predicate<SAMRecord> filterSamRecords=new Predicate<SAMRecord>()
				{
				@Override
				public boolean apply(SAMRecord rec)
					{
					if(rec.getReadUnmappedFlag()) return false;
					if(rec.isSecondaryOrSupplementary()) return false;
					if(rec.getDuplicateReadFlag()) return false;
					if(rec.getReadFailsVendorQualityCheckFlag()) return false;
					Cigar cigar=rec.getCigar();
					if(cigar==null ||  cigar.numCigarElements()<2) return false;
					
					boolean found_S=false;
					for(int side=0;side<2;++side)
						{
						CigarElement ce=cigar.getCigarElement(side==0?0:cigar.numCigarElements()-1);
						//read must be clipped on 5' or 3' with a good length
						if(!ce.getOperator().equals(CigarOperator.S))  continue;
						found_S=true;
						break;
						}
					if(!found_S) return false;
					SAMReadGroupRecord g=rec.getReadGroup();
					if(g==null || g.getSample()==null || g.getSample().isEmpty()) return false;		
					return true;
					}
				};
			final FilteringIterator<SAMRecord> forwardIterator=new FilteringIterator<SAMRecord>(merginIter, filterSamRecords);
			for(;;)
				{
				
				SAMRecord rec=null;
				if(forwardIterator.hasNext())
					{
					rec=forwardIterator.next();
					progress.watch(rec);
					if(intervals!=null && !intervals.containsOverlapping(rec.getReferenceName(),rec.getAlignmentStart(),rec.getAlignmentEnd())) continue;
					}
				//need to flush buffer ?
				if( rec==null ||
					(!buffer.isEmpty() && !buffer.getLast().getReferenceIndex().equals(rec.getReferenceIndex())) ||
					(!buffer.isEmpty() && buffer.getLast().getUnclippedEnd()+readLength < rec.getUnclippedStart())
					)
					{
					if(!buffer.isEmpty())
						{
						int chromStart=buffer.getFirst().getUnclippedStart();
						int chromEnd=buffer.getFirst().getUnclippedEnd();
						for(SAMRecord sam:buffer)
							{
							chromStart=Math.min(chromStart, sam.getUnclippedStart());
							chromEnd=Math.max(chromEnd, sam.getUnclippedEnd());
							}
						final int winShift=5;
						for(int pos=chromStart;pos+winShift<=chromEnd;pos+=winShift)
							{
							int count_big_clip[]=new int[]{0,0};
							//int max_depth[]=new int[]{0,0};
							List<Genotype> genotypes=new ArrayList<Genotype>();	
							Set<Allele> all_alleles=new HashSet<Allele>();
							all_alleles.add(reference_allele);
							boolean found_one_depth_ok=false;
							int sum_depth=0;
							int samples_with_high_depth=0;
							for(String sample: sample2input.keySet())
								{
								GenotypeBuilder gb=new GenotypeBuilder(sample);
								int count_clipped[]=new int[]{0,0};
								Set<Allele> sample_alleles=new HashSet<Allele>(3);
								for(int side=0;side<2;++side)
									{
									for(SAMRecord sam:buffer)
										{
										if(!sam.getReadGroup().getSample().equals(sample)) continue;
										Cigar cigar=sam.getCigar();
										CigarElement ce=cigar.getCigarElement(side==0?0:cigar.numCigarElements()-1);
										if(!ce.getOperator().equals(CigarOperator.S))  continue;
										int clipStart=(side==0?sam.getUnclippedStart():sam.getAlignmentEnd()+1);
										int clipEnd=(side==0?sam.getAlignmentStart()-1:sam.getUnclippedEnd());
										if((pos+winShift< clipStart || pos> clipEnd)) continue;
										count_clipped[side]++;
										if(ce.getLength()>=this.min_clip_length)
											{
											count_big_clip[side]++;
											}
										sample_alleles.add(alternate_alleles[side]);
										gb.attribute("CN"+(side==0?5:3),count_clipped[side]);
										}
									}
								//if(!(found_one_big_clip[0] || found_one_big_clip[1])) continue;
								if(count_clipped[0]+count_clipped[1]==0) continue;
								if((count_clipped[0]+count_clipped[1]) > min_depth)
									{
									found_one_depth_ok=true;
									++samples_with_high_depth;
									}
								sum_depth+=(count_clipped[0]+count_clipped[1]);
								gb.alleles(new ArrayList<Allele>(sample_alleles));
								all_alleles.addAll(sample_alleles);
								gb.DP(count_clipped[0]+count_clipped[1]);
								genotypes.add(gb.make());
								}
							
							if(all_alleles.size()==1)
								{
								continue;//all homozygotes
								}
							if(!found_one_depth_ok)
								{
								continue;
								}
							if(!(count_big_clip[0]>=1 || count_big_clip[1]>=1) )
								{
								continue;
								}
							Map<String,Object> atts=new HashMap<String,Object>();
							atts.put("COUNT_SAMPLES",samples_with_high_depth);
							atts.put(VCFConstants.DEPTH_KEY, sum_depth);
							VariantContextBuilder vcb=new VariantContextBuilder();
							vcb.chr(buffer.getFirst().getReferenceName());
							vcb.start(pos);
							vcb.stop(pos+winShift);
							vcb.alleles(all_alleles);
							vcb.attributes(atts);
							vcb.genotypes(genotypes);
							w.add(vcb.make());
							}
						
						buffer.clear();
						}
					if(rec==null)
						{
						break;
						}
					}
				
				buffer.add(rec);
				
				
				

				}
			
			merginIter.close();
			progress.finish();
			return 0;
			}

		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			for(Input input:inputs)
				{
				CloserUtil.close(input);
				}
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new SamFindClippedRegions().instanceMainWithExit(args);
		}
	}
