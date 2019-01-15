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
package com.github.lindenb.jvarkit.tools.structvar;

import java.io.BufferedReader;
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

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFHeaderVersion;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFSimpleHeaderLine;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.MergingSamRecordIterator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamFileHeaderMerger;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.SequenceUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;

import java.util.function.Predicate;

import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
/**
BEGIN_DOC


END_DOC

 */
@Program(name="samfindclippedregions")
public class SamFindClippedRegions extends Launcher
	{
	private static final Logger LOG=Logger.build(SamFindClippedRegions.class).make();

	@Parameter(names="-B",description="bed file")
	private File bedFile = null;

	@Parameter(names="-c",description="min size of clipped read")
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
				if(filter.test(o))
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
		SamReader samFileReaderScan;
		SAMFileHeader header;
		String sampleName=null;
		Input(File bamFile)
			{
			LOG.info("Reading from "+bamFile);
			this.bamFile=bamFile;
			
			SamReaderFactory srf=SamReaderFactory.make().
					validationStringency(ValidationStringency.LENIENT);

			
			this.samFileReaderScan=srf.open(this.bamFile);
			this.header=this.samFileReaderScan.getFileHeader();
			for(SAMReadGroupRecord g:this.header.getReadGroups())
				{
				if(g.getSample()==null || g.getSample().isEmpty())
					{
					LOG.warning("Read group "+g.getId()+" has no sample in "+bamFile);
					continue;
					}
				if(sampleName!=null && !sampleName.equals(g.getSample()))
					{	
					throw new RuntimeException("Sorry one bam can only contain one sample and I got "+sampleName+" and "+g.getSample()+" in "+bamFile);
					}
				sampleName=g.getSample();
				}
			if(sampleName==null)
				{
				throw new RuntimeException("Bam "+this.bamFile+" doesn't contain a sample defined in header/read-group");
				}
			LOG.info("In "+bamFile+" Sample is "+sampleName);
			}
		
		@Override
		public void close() throws IOException {
			CloserUtil.close(samFileReaderScan);
			}
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
	public int doWork(List<String> args) {
		int readLength=150;
		
		if(args.isEmpty())
			{
			LOG.error("illegal.number.of.arguments");
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
			for(final String filename:args)
				{
				Input input=new Input(new File(filename));
				//input.index=inputs.size();
				inputs.add(input);
				
				if(sample2input.containsKey(input.sampleName))
					{
					LOG.error("Duplicate sample "+input.sampleName+" in "+input.bamFile+" and "+sample2input.get(input.sampleName).bamFile);
					return -1;
					}
				sample2input.put(input.sampleName, input);
				if(dict==null)
					{
					dict=input.header.getSequenceDictionary();
					}
				else if(!SequenceUtil.areSequenceDictionariesEqual(dict, input.header.getSequenceDictionary()))
					{
					LOG.error("Found more than one dictint sequence dictionary");
					return -1;
					}
					
				}
			LOG.info("Sample N= "+sample2input.size());
			
			
			/* create merged iterator */
			List<SAMFileHeader> headers=new ArrayList<SAMFileHeader>(sample2input.size());
			for(Input input:inputs) headers.add(input.header);
			SamFileHeaderMerger headerMerger=new SamFileHeaderMerger(SortOrder.coordinate,headers, true);
			
			List<SamReader> readers=new ArrayList<SamReader>(sample2input.size());
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
			

			
			final IntervalTreeMap<Boolean> intervals= new IntervalTreeMap<>();
			//w=swf.make(header, System.out);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
			
			if(bedFile!=null)
				{
				final BedLineCodec bedLineCodec=new BedLineCodec();
				LOG.info("Reading "+bedFile);
				BufferedReader r=IOUtils.openFileForBufferedReading(bedFile);
				String line;
				while((line=r.readLine())!=null)
					{
					BedLine bedLine = bedLineCodec.decode(line);
					if(bedLine==null) continue;

					if(dict!=null && dict.getSequence(bedLine.getContig())==null)
						{
						LOG.warning("undefined chromosome  in "+bedFile+" "+line);
						continue;
						}
					intervals.put(bedLine.toInterval(), true);
					}
				CloserUtil.close(r);
				}
			
			
			
			LinkedList<SAMRecord> buffer=new LinkedList<SAMRecord>();
			
			final Predicate<SAMRecord> filterSamRecords=new Predicate<SAMRecord>()
				{
				@Override
				public boolean test(SAMRecord rec)
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
					if(intervals!=null && !intervals.containsOverlapping(
							new Interval(rec.getReferenceName(),rec.getAlignmentStart(),rec.getAlignmentEnd()))) continue;
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
			LOG.error(err);
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
