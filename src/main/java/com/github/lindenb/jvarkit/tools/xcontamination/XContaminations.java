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
import htsjdk.samtools.filter.SamRecordFilter;
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

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamFilterParser;
import com.github.lindenb.jvarkit.util.illumina.ShortReadName;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

/**
BEGIN_DOC

## Example

```bash
$ find . -type f -name "*.bam" > bam.list
$  head -n 10000 variant.vcf | java -jar dist/xcontaminations.jar - bam.list > out.tsv
$ verticalize out.tsv


>>> 2
$1       #Machine:FlowCell:Run:Lane-1 : HISEQ10:C3FBPACXX:0:4
$2                            sample1 : B00G5V9
$3        Machine:FlowCell:Run:Lane-2 : HISEQ10:C486PACXX:0:3
$4                            sample2 : B00G7LK
$5                          same.lane : 0
$6   reads_sample1_supporting_sample1 : 26392
$7   reads_sample1_supporting_sample2 : 70
$8    reads_sample1_supporting_others : 40
$9   reads_sample2_supporting_sample2 : 21473
$10  reads_sample2_supporting_sample1 : 39
$11    reads_sample2_supporting_other : 31
<<< 2

(...)

>>> 9
$1       #Machine:FlowCell:Run:Lane-1 : HISEQ5:C3FV0ACXX:0:7
$2                            sample1 : B00G738
$3        Machine:FlowCell:Run:Lane-2 : HISEQ5:C3FV0ACXX:0:7
$4                            sample2 : B00G754
$5                          same.lane : 1
$6   reads_sample1_supporting_sample1 : 10209
$7   reads_sample1_supporting_sample2 : 23
$8    reads_sample1_supporting_others : 15
$9   reads_sample2_supporting_sample2 : 9054
$10  reads_sample2_supporting_sample1 : 32
$11    reads_sample2_supporting_other : 9
<<< 9

```

generating in parallel:

```make
CHROMS=1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22

.PHONY:all

define xcont

$$(addprefix tmp.xcont.,$$(addsuffix .tsv.gz,$(1))) :
        bcftools view -r "$(1)" unifiedgenotyper.vcf.gz -Tcapture.bed | java -Xmx1g -jar xcontaminations.jar - bal.list  | gzip --best > $$(addsuffix .tmp.gz,$$@) && mv  
$$(addsuffix .tmp.gz,$$@) $$@


endef

all: $(foreach C,${CHROMS},$(addprefix tmp.xcont.,$(addsuffix .tsv.gz,${C})))
        $(foreach I,0 1, gunzip -c $^ |  awk -F '       ' '($$5==$I)'  |awk -F '        ' 'BEGIN {T=0;N=0;} {for(i=6;i<=NF;++i) T+=int($$i); N+=int($$7); N+=int($$10);} E
ND { printf("%f\n",N/T);}'; )

$(foreach C,${CHROMS},$(eval $(call xcont,$C)))
```



END_DOC

 */

@Program(name="xcontaminations",description="For @AdrienLeger2 : cross contamination between samples in same lane")
public class XContaminations extends Launcher
	{
	private static final Logger LOG=Logger.build(XContaminations.class).make();
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-filter","--filter"},description=SamFilterParser.FILTER_DESCRIPTION,converter=SamFilterParser.StringConverter.class)
	private SamRecordFilter filter  = SamFilterParser.buildDefault();

	
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
	
	
	
	public void addBamFile(File bamFile)
		{
		this.bamFiles.add(bamFile);
		}
	
	@Override
	public int doWork(List<String> args) {
		
		if(args.size()<2)
			{
			LOG.error("Illegal Number of args");
			return -1;
			}
		
		for(String f:IOUtils.unrollFiles(args.subList(1, args.size())))
			{
			addBamFile(new File(f));
			}
			
		if(this.bamFiles.isEmpty())
			{
			LOG.error("Undefined BAM file(s)");
			return -1;
			}	
		
		SAMRecordIterator iter=null;
		VcfIterator in=null;
		Map<String,SamReader> sample2samReader=new HashMap<>();

		try {
			SamReaderFactory srf=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
			
			if(args.get(0).equals("-"))
				{
				in = super.openVcfIterator(null);
				}
			else
				{
				in = super.openVcfIterator(args.get(0));
				}
			
			
			VCFHeader vcfHeader=in.getHeader();
			SAMSequenceDictionary dict1=vcfHeader.getSequenceDictionary();
			if(dict1==null)
				{
				LOG.error("VCF is missing a SAM sequence dictionary");
				return -1;
				}
			
			Set<String> sampleNames= new HashSet<>(vcfHeader.getSampleNamesInOrder());
			if( sampleNames.isEmpty())
				{
				LOG.error("VCF contains no sample");
				return -1;
				}
			
			for(File bamFile:this.bamFiles)
				{
				LOG.info("Opening "+bamFile);
				SamReader samReader=srf.open(bamFile);
				SAMFileHeader samHeader= samReader.getFileHeader();
				SAMSequenceDictionary dict2=samHeader.getSequenceDictionary();
				if(dict2==null)
					{
					samReader.close();
					LOG.error("BAM is missing a SAM sequence dictionary");
					return -1;
					}
				
				if(!SequenceUtil.areSequenceDictionariesEqual(dict1, dict2))
					{
					samReader.close();
					LOG.error("VCF/BAM Not the same dictionaries");
					return -1;
					}
				
				if(!samReader.hasIndex())
					{
					samReader.close();
					LOG.error("sam is not index");
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
						LOG.error("Cannot handle more than one sample/bam");
						return -1;
						}
					}
				if(sampleName==null)
					{
					samReader.close();
					LOG.error("No sample in "+bamFile);
					continue;//skip this bam
					}
				if(!sampleNames.contains(sampleName))
					{
					samReader.close();
					LOG.error("Not in VCF header: sample "+sampleName+" "+bamFile);
					continue;//skip this bam
					}
				if(sample2samReader.containsKey(sampleName))
					{
					samReader.close();
					LOG.error("Cannot handle more than one bam/sample");
					return -1;
					}
				
				sample2samReader.put(sampleName, samReader);
				}
			
			if(sample2samReader.size()<2)
				{
				LOG.error("Not engough BAM/samples. Expected at least two valid BAMs");
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
							ctx.getContig(),
							ctx.getStart(),
							ctx.getEnd(),
							false
							);
					
					while(iter.hasNext())
						{
						SAMRecord record= iter.next();
						if(record.getReadUnmappedFlag()) continue;
						if(this.filter.filterOut(record)) continue;
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
							LOG.info("No a valid read name "+record.getReadName());
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
									if( contaminationTable.size()%10000==0) LOG.info("n(pairs)=" + contaminationTable.size() ); 
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
			PrintWriter pw= super.openFileOrStdoutAsPrintWriter(this.outputFile);
			
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
			LOG.error(e);
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

	
	
	public static void main(String[] args)
		{
		new XContaminations().instanceMainWithExit(args);
		}
	}
