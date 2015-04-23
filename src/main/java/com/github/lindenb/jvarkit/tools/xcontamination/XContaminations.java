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
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.xcontamination;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

import java.io.File;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.knime.AbstractKnimeApplication;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.illumina.ShortReadName;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class XContaminations extends AbstractKnimeApplication
	{
	private Set<File> bamFiles=new HashSet<File>();
	
	private static class SampleAlleles
		{
		long reads_sample1_supporting_sample1=0;
		long reads_sample1_supporting_sample2=0;
		long reads_sample1_supporting_other=0;
		long reads_sample2_supporting_sample1=0;
		long reads_sample2_supporting_sample2=0;
		long reads_sample2_supporting_other=0;
		}
	
	private static class SamplePair
		{
		SequencerFlowCellRunLaneSample sample1;
		SequencerFlowCellRunLaneSample sample2;
		SamplePair(SequencerFlowCellRunLaneSample s1,SequencerFlowCellRunLaneSample s2)
			{
			if(s1.sampleName.compareTo(s2.sampleName)<0)
				{
				sample1=s1;
				sample2=s2;
				}
			else
				{
				sample1=s2;
				sample2=s1;
				}
			}
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + sample1.hashCode();
			result = prime * result + sample2.hashCode();
			return result;
			}

		@Override
		public boolean equals(Object obj) {
			if (this == obj) return true;
			if (obj == null) return false;
			if (getClass() != obj.getClass()) return false;
			SamplePair other = (SamplePair) obj;
			if (!sample1.equals(other.sample1)) return false;
			if (!sample2.equals(other.sample2)) return false;
			return true;
			}
		@Override
		public String toString() {
			return "("+sample1+"/"+sample2+")";
			}
		}
	
	private static class SequencerFlowCellRunLaneSample
		{
		String machine;
		String flowCell;
		int run;
		int lane;
		String sampleName;
		
		SequencerFlowCellRunLaneSample(ShortReadName name,String sampleName)
			{
			this.machine=name.getInstrumentName();
			this.flowCell=name.getFlowCellId();
			this.run=Math.max(name.getRunId(),0);
			this.lane=name.getFlowCellLane();
			this.sampleName=sampleName;
			}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + flowCell.hashCode();
			result = prime * result + lane;
			result = prime * result + machine.hashCode();
			result = prime * result + sampleName.hashCode();
			result = prime * result + run;
			return result;
			}

		@Override
		public boolean equals(Object obj) {
			if (this == obj) return true;
			if (obj == null) return false;
			if (getClass() != obj.getClass()) return false;
			SequencerFlowCellRunLaneSample other = (SequencerFlowCellRunLaneSample) obj;
			if (lane != other.lane) return false;
			if (run != other.run) return false;
			if (!flowCell.equals(other.flowCell)) return false;
			if (!machine.equals(other.machine)) return false;
			if (!sampleName.equals(other.sampleName)) return false;
			return true;
			}
		
		public String getSequencingLabel()
			{
			return machine+":"+flowCell+":"+run+":"+lane;
			}
		
		@Override
		public String toString() {
			return getSequencingLabel()+":"+sampleName;
			}
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return DEFAULT_WIKI_PREFIX+"XContaminations";
		}
	@Override
	public String getProgramDescription() {
		return "For @AdrienLeger2 : cross contamination between samples in same lane";
		}
	
	public void addBamFile(File bamFile)
		{
		this.bamFiles.add(bamFile);
		}
	
	@Override
	public int executeKnime(List<String> args)
		{
		if(args.size()<2)
			{
			error("Illegal Number of args");
			return -1;
			}
		
		for(String f:IOUtils.unrollFiles(args.subList(1, args.size())))
			{
			addBamFile(new File(f));
			}
			
		if(this.bamFiles.isEmpty())
			{
			error("Undefined BAM file(s)");
			return -1;
			}	
		
		SAMRecordIterator iter=null;
		VcfIterator in=null;
		Map<String,SamReader> sample2samReader=new HashMap<>();

		try {
			SamReaderFactory srf=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
			
			if(args.get(0).equals("-"))
				{
				in = VCFUtils.createVcfIteratorStdin();
				}
			else
				{
				in = VCFUtils.createVcfIteratorFromFile(new File(args.get(0)));
				}
			
			
			VCFHeader vcfHeader=in.getHeader();
			SAMSequenceDictionary dict1=vcfHeader.getSequenceDictionary();
			if(dict1==null)
				{
				error("VCF is missing a SAM sequence dictionary");
				return -1;
				}
			
			Set<String> sampleNames= new HashSet<>(vcfHeader.getSampleNamesInOrder());
			if( sampleNames.isEmpty())
				{
				error("VCF contains no sample");
				return -1;
				}
			
			for(File bamFile:this.bamFiles)
				{
				info("Opening "+bamFile);
				SamReader samReader=srf.open(bamFile);
				SAMFileHeader samHeader= samReader.getFileHeader();
				SAMSequenceDictionary dict2=samHeader.getSequenceDictionary();
				if(dict2==null)
					{
					samReader.close();
					error("BAM is missing a SAM sequence dictionary");
					return -1;
					}
				
				if(!SequenceUtil.areSequenceDictionariesEqual(dict1, dict2))
					{
					samReader.close();
					error("VCF/BAM Not the same dictionaries");
					return -1;
					}
				
				if(!samReader.hasIndex())
					{
					samReader.close();
					error("sam is not index");
					return -1;
					}
				String sampleName=null;
				for(SAMReadGroupRecord rgr:samHeader.getReadGroups())
					{
					String s=rgr.getSample();
					if(s==null ) continue;
					if(sampleName==null)
						{
						sampleName=s;
						}
					else if(!sampleName.equals(s))
						{
						samReader.close();
						error("Cannot handle more than one sample/bam");
						return -1;
						}
					}
				if(sampleName==null)
					{
					samReader.close();
					error("No sample in "+bamFile);
					continue;//skip this bam
					}
				if(!sampleNames.contains(sampleName))
					{
					samReader.close();
					error("Not in VCF header: sample "+sampleName+" "+bamFile);
					continue;//skip this bam
					}
				if(sample2samReader.containsKey(sampleName))
					{
					samReader.close();
					error("Cannot handle more than one bam/sample");
					return -1;
					}
				
				sample2samReader.put(sampleName, samReader);
				}
			
			if(sample2samReader.size()<2)
				{
				error("Not engough BAM/samples. Expected at least two valid BAMs");
				return -1;
				}
			
			sampleNames.retainAll(sample2samReader.keySet());
			
			Map<SamplePair,SampleAlleles> contaminationTable=new HashMap<>();
			
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict1);
			while(in.hasNext())
				{
				VariantContext ctx= progress.watch(in.next());
				if(!ctx.isSNP() || !ctx.isBiallelic() || ctx.isSymbolic()) continue;
				
				
				boolean isWorthScanning=false;
				/* scan Reads for this Sample */
				for(String s1: sampleNames)
					{
					Genotype g1=ctx.getGenotype(s1);
					if(g1==null || !g1.isHom()) continue;
					/* scan Reads for this Sample */
					for(String s2: sampleNames)
						{
						if(s1.compareTo(s2)>=0) continue;
						Genotype g2=ctx.getGenotype(s2);
						if(g2==null || !g2.isHom()) continue;
						if(g1.sameGenotype(g2)) continue;
						isWorthScanning=true;
						break;
						}
					if(isWorthScanning) break;
					}
				
				if(!isWorthScanning) continue;

				//Set<SequencerFlowCellRunLane> sequencerFlowCellRunLaneInThisContext=new HashSet<>();
				Map<SequencerFlowCellRunLaneSample,Counter<Allele>> sample2allelesCount=new HashMap<>();
				
				/* scan Reads for this Sample */
				for(String sampleName: sampleNames)
					{
					//sample name is not in vcf header
					if(!sample2samReader.containsKey(sampleName))
						continue;
					

					
					SamReader samReader = sample2samReader.get(sampleName);
					if(samReader==null) continue;
					
					Genotype genotype = ctx.getGenotype(sampleName);
					if(genotype==null || !genotype.isHom()) continue;
					iter = samReader.query(
							ctx.getChr(),
							ctx.getStart(),
							ctx.getEnd(),
							false
							);
					
					while(iter.hasNext())
						{
						SAMRecord record= iter.next();
						if(record.getReadUnmappedFlag()) continue;
						if(record.isSecondaryOrSupplementary()) continue;
						if(record.getDuplicateReadFlag()) continue;
						if(record.getMappingQuality()==0 || record.getMappingQuality()==255) continue;
						if(record.getReadPairedFlag())
							{
							if(!record.getProperPairFlag()) continue;
							}
						SAMReadGroupRecord srgr = record.getReadGroup();
						
						//not current sample
						if(srgr==null) continue;
						if(!sampleName.equals(srgr.getSample())) continue;
						
						ShortReadName readName = ShortReadName.parse(record);
						if(!readName.isValid())
							{
							info("No a valid read name "+record.getReadName());
							continue;
							}
						
						
						SequencerFlowCellRunLaneSample sequencerFlowCellRunLaneSample=
								new SequencerFlowCellRunLaneSample(readName, sampleName);
						
						
						
						Cigar cigar=record.getCigar();
						if(cigar==null) continue;
						byte readSeq[]=record.getReadBases();
						if(readSeq==null) continue;
						int refPos1 = record.getAlignmentStart();
						int readPos = 0; 
						for(CigarElement ce: cigar.getCigarElements())
							{
							//beyond variant position ?
							if(refPos1>ctx.getStart()) break;
							CigarOperator op=ce.getOperator();
							switch(op)
								{
								case I: readPos+=ce.getLength(); break;
								case N://threw
								case D: refPos1+=ce.getLength(); break;
								case S: readPos+=ce.getLength(); break;
								case H: break;
								case P: break;
								case M: case EQ: case X:
									{
									for(int i=0;i< ce.getLength();++i)
										{
										if( refPos1 == ctx.getStart())
											{	
											byte base =  readSeq[readPos];
											if(!(base=='N' || base=='n'))
												{
												Counter<Allele> sampleAlleles= sample2allelesCount.get(sequencerFlowCellRunLaneSample);
												if(sampleAlleles==null)
													{
													sampleAlleles=new Counter<Allele>();
													sample2allelesCount.put(sequencerFlowCellRunLaneSample, sampleAlleles);

													}

												sampleAlleles.incr(
													Allele.create(base,false)
													);

												}
											}
										refPos1++;
										readPos++;
										}
									break;
									}
								default: throw new IllegalStateException();
								}
							}
						}
					iter.close();
					iter=null;
					}/* end scan reads for this sample */
				/* sum-up data for this SNP */
			
				for(String sample1: sampleNames)
					{
					Genotype g1= ctx.getGenotype(sample1);
					if(g1==null || !g1.isHom()) continue;
					Allele a1 = Allele.create(g1.getAllele(0).getBases(),false);
					
					
					for(String sample2: sampleNames)
						{
						if(sample1.compareTo(sample2)>=0) continue;
						Genotype g2= ctx.getGenotype(sample2);
						if(g2==null || !g2.isHom() || g2.sameGenotype(g1)) continue;
						Allele a2 = Allele.create(g2.getAllele(0).getBases(),false);
						
						for(SequencerFlowCellRunLaneSample sfcr1:sample2allelesCount.keySet())
							{
							if(!sfcr1.sampleName.equals(sample1)) continue;
							for(SequencerFlowCellRunLaneSample sfcr2:sample2allelesCount.keySet())
								{
								if(!sfcr2.sampleName.equals(sample2)) continue;
								
								SamplePair samplePair=new SamplePair(sfcr1, sfcr2);
								
								Counter<Allele> counter1 =  sample2allelesCount.get(sfcr1);
								if(counter1==null) continue;
								Counter<Allele> counter2 =  sample2allelesCount.get(sfcr2);
								if(counter2==null) continue;
								
								
								SampleAlleles sampleAlleles= contaminationTable.get(samplePair);
								if(sampleAlleles==null)
									{
									sampleAlleles=new SampleAlleles();
									contaminationTable.put(samplePair,sampleAlleles);
									if( contaminationTable.size()%10000==0) info("n(pairs)=" + contaminationTable.size() ); 
									}
								
								for(Allele allele: counter1.keySet())
									{
									long n = counter1.count(allele);
									if(allele.equals(a1))
										{
										sampleAlleles.reads_sample1_supporting_sample1 += n;
										}
									else if(allele.equals(a2))
										{
										sampleAlleles.reads_sample1_supporting_sample2 += n;
										}
									else
										{
										sampleAlleles.reads_sample1_supporting_other += n;
										}
									}
								
								for(Allele allele: counter2.keySet())
									{
									long n = counter2.count(allele);
									if(allele.equals(a2))
										{
										sampleAlleles.reads_sample2_supporting_sample2 += n;
										}
									else if(allele.equals(a1))
										{
										sampleAlleles.reads_sample2_supporting_sample1 += n;
										}
									else
										{
										sampleAlleles.reads_sample2_supporting_other += n;
										}
									}
								
								}
							
							
							}
						}
				
					}
				
				}
			progress.finish();
			boolean somethingPrinted=false;
			PrintWriter pw= null;
			if(getOutputFile()==null)
				{
				pw=new PrintWriter(System.out);
				}
			else
				{
				pw=new PrintWriter(getOutputFile());
				}
			/* we're done, print the result */
			pw.print("#");
			pw.print("Machine:FlowCell:Run:Lane-1\tsample1");
			pw.print("\tMachine:FlowCell:Run:Lane-2\tsample2");
			pw.print("\tsame.lane");
			pw.print("\treads_sample1_supporting_sample1");
			pw.print('\t');
			pw.print("reads_sample1_supporting_sample2");
			pw.print('\t');
			pw.print("reads_sample1_supporting_others");
			pw.print('\t');
			pw.print("reads_sample2_supporting_sample2");
			pw.print('\t');
			pw.print("reads_sample2_supporting_sample1");
			pw.print('\t');
			pw.print("reads_sample2_supporting_other");

			pw.println();
			for(SamplePair pair : contaminationTable.keySet())
				{
				SampleAlleles sampleAlleles = contaminationTable.get(pair);
				if(sampleAlleles==null) continue;
				
				pw.print(pair.sample1.getSequencingLabel());
				pw.print('\t');
				pw.print(pair.sample1.sampleName);
				pw.print('\t');
				pw.print(pair.sample2.getSequencingLabel());
				pw.print('\t');
				pw.print(pair.sample2.sampleName);
				pw.print('\t');
				pw.print(pair.sample1.getSequencingLabel().equals(pair.sample2.getSequencingLabel())?1:0);
				pw.print('\t');
				pw.print(sampleAlleles.reads_sample1_supporting_sample1);
				pw.print('\t');
				pw.print(sampleAlleles.reads_sample1_supporting_sample2);
				pw.print('\t');
				pw.print(sampleAlleles.reads_sample1_supporting_other);
				pw.print('\t');
				pw.print(sampleAlleles.reads_sample2_supporting_sample2);
				pw.print('\t');
				pw.print(sampleAlleles.reads_sample2_supporting_sample1);
				pw.print('\t');
				pw.print(sampleAlleles.reads_sample2_supporting_other);
				pw.println();
				somethingPrinted=true;
						
					
				
				
				
				}
			pw.flush();
			pw.close();
			
			if(!somethingPrinted)
				{
				System.err.println("Warning: NO output");
				}
			
			return 0;
			}
		catch (Exception e) {
			error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			CloserUtil.close(iter);
			for(SamReader samReader:sample2samReader.values())
				CloserUtil.close(samReader);
			sample2samReader.clear();
			}
		
		}

	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "o:"))!=-1)
			{
			switch(c)
				{
				case 'o': setOutputFile(opt.getOptArg()); break;
				default: 
					{
					switch(handleOtherOptions(c, opt, args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default: break;
						}
					}
				}
			}
		return mainWork(opt.getOptInd(), args);
		}
	
	public static void main(String[] args)
		{
		new XContaminations().instanceMainWithExit(args);
		}
	}
