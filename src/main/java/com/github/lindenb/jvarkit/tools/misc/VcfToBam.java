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
package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeSet;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import htsjdk.variant.vcf.VCFIterator;

/*
 
BEGIN_DOC

##Example

```bash
$ java -jar dist-1.139/vcf2bam.jar -R ref.fa  samtools.vcf.gz 2> /dev/null | grep -v "100="



@HD	VN:1.5	SO:unsorted
@SQ	SN:seg1	LN:5101
@SQ	SN:seg2	LN:4000
@RG	ID:NA12878	SM:NA12878	LB:illumina	DS:NA12878
@RG	ID:NA12891	SM:NA12891	LB:illumina	DS:NA12891
@RG	ID:NA12892	SM:NA12892	LB:illumina	DS:NA12892
@CO	Generated with -R /home/lindenb/src/ngsxml/test/ref/ref.fa /home/lindenb/src/ngsxml/OUT/Projects/Proj1/VCF/samtools/Proj1.samtools.vcf.gz
0000003609	83	seg1	1449	60	99=1X	=	949	-500	TGCTTTAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12891	NM:i:1
0000003610	147	seg1	1449	60	99=1X	=	949	-500	TGCTTTAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12891	NM:i:1
0000003611	147	seg1	1449	60	99=1X	=	949	-500	TGCTTTAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12892	NM:i:1
0000003612	83	seg1	1449	60	99=1X	=	949	-500	TGCTTTAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12892	NM:i:1
0000003615	83	seg1	1450	60	98=1X1=	=	950	-500	GCTTTAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTCT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12891	NM:i:1
0000003616	147	seg1	1450	60	98=1X1=	=	950	-500	GCTTTAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTCT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12891	NM:i:1
0000003617	147	seg1	1450	60	98=1X1=	=	950	-500	GCTTTAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTCT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12892	NM:i:1
0000003627	147	seg1	1452	60	96=1X3=	=	952	-500	TTTAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTCTCT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12891	NM:i:1
0000003629	147	seg1	1452	60	96=1X3=	=	952	-500	TTTAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTCTCT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12892	NM:i:1
0000003630	147	seg1	1452	60	96=1X3=	=	952	-500	TTTAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTCTCT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12892	NM:i:1
0000003633	147	seg1	1453	60	95=1X4=	=	953	-500	TTAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTCTCTC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12891	NM:i:1
0000003635	83	seg1	1453	60	95=1X4=	=	953	-500	TTAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTCTCTC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12892	NM:i:1
0000003640	147	seg1	1454	60	94=1X5=	=	954	-500	TAGAAAGCATTCCAAAATCTCTTACCAGTTTTATCTCCTATGAAAGTCCTTCACACTTTCTCTCATTTAAACTTTATTGCATTTTCCTCACTTTCTCTCA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	RG:Z:NA12891	NM:i:1

```


END_DOC

*/

@Program(name="vcf2bam",description="vcf to bam",keywords={"ref","vcf","bam"})
public class VcfToBam extends Launcher
	{
	private static final Logger LOG=Logger.build(VcfToBam.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile=null;
	
	@Parameter(names={"-r","-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private File faidx=null;

	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();
	
	private IndexedFastaSequenceFile indexedFastaSequenceFile;
	@Parameter(names="--fragmentsize",description="fragment size")
	private int fragmentSize=600;
	@Parameter(names="--readsize",description="read size")
	private int readSize=100;
	
	
	
	
	
	private void run(VCFIterator vcfIterator) throws IOException
		{
		 long id_generator=0L;
		 SAMFileWriter samFileWriter =null;
		
		final VCFHeader header = vcfIterator.getHeader();
		SAMSequenceDictionary dict = header.getSequenceDictionary();
		if(dict==null) throw new IOException("Sequence Dictionary missing in VCF");
		if(!SequenceUtil.areSequenceDictionariesEqual(dict, indexedFastaSequenceFile.getSequenceDictionary()))
			{
			LOG.warning(
					"REF/VCF: not the same Sequence dictonary "+dict+" "+ indexedFastaSequenceFile.getSequenceDictionary()
					);
			}
		final SAMFileHeader samHeader=new SAMFileHeader();
		samHeader.setSequenceDictionary(dict);
		for(final String sample:new TreeSet<>(header.getSampleNamesInOrder()))
			{
			final SAMReadGroupRecord rg= new SAMReadGroupRecord(sample);
			rg.setSample(sample);
			rg.setLibrary(sample);
			rg.setDescription(sample);
			rg.setLibrary("illumina");
			samHeader.addReadGroup(rg);
			}
		samHeader.addComment("Generated with "+getProgramCommandLine());
		samHeader.setSortOrder(SortOrder.unsorted);
		samFileWriter= this.writingBamArgs.
				setReferenceFile(this.faidx).
				openSAMFileWriter(this.outputFile, samHeader, true);
		
		/* looping over sequences */
		for(final SAMSequenceRecord ssr: dict.getSequences())
			{
			final LinkedList<VariantContext> variantBuffer = new LinkedList<>();
			GenomicSequence genomicSequence=new GenomicSequence(this.indexedFastaSequenceFile, ssr.getSequenceName());
			int x=1;
			while(x+this.fragmentSize <= ssr.getSequenceLength())
				{
				// pop on left
				while(!variantBuffer.isEmpty() && (variantBuffer.getFirst().getStart()+this.fragmentSize*2) < x)
					{
					variantBuffer.removeFirst();
					}
				//refill buffer
				while(vcfIterator.hasNext() )
					{
					//buffer is already by enough
					if(!variantBuffer.isEmpty() && variantBuffer.getLast().getStart()> x+2*fragmentSize)
						{
						break;
						}
					
					final VariantContext ctx = vcfIterator.peek();
					if(ctx==null) break;
					if(ctx.isIndel() || !ctx.isVariant())
						{
						LOG.warning("Cannot handle "+ctx.getReference()+"/"+ctx.getAlternateAlleles());
						vcfIterator.next();//consumme
						continue;
						}
					
					final SAMSequenceRecord variantContig = dict.getSequence( ctx.getContig());
					if(variantContig == null) throw new IOException("Unknown contig "+ctx.getContig());
					if(variantContig.getSequenceIndex()< ssr.getSequenceIndex())
						{
						//just consumme !
						// https://github.com/lindenb/jvarkit/issues/86#issuecomment-326986654
						vcfIterator.next();//consumme
						continue;
						}	
					else if(variantContig.getSequenceIndex() > ssr.getSequenceIndex())
						{
						break;
						}
					else
						{
						variantBuffer.add(vcfIterator.next());
						}
					}
				
				if(variantBuffer.isEmpty())
					{
					LOG.info("no variant ?");
					//no variant on this chromosome
					break;
					}
				if(x  < variantBuffer.getFirst().getStart()-2*fragmentSize)
					{
					x= variantBuffer.getFirst().getStart()-2*fragmentSize;
					}
				
				for(int depth =0;depth<1;++depth)
					{
					for(String sample:header.getSampleNamesInOrder())
						{
						//loop over genomic strand
						for(int g_strand=0;g_strand<2;++g_strand)
							{
							
							SAMRecord records[]={
									new SAMRecord(samHeader),	
									new SAMRecord(samHeader)
									};
							String readName=String.format("%010d",++id_generator);
							for(int R=0;R<2;++R)
								{
								records[R].setReferenceName(ssr.getSequenceName());
								records[R].setReadPairedFlag(true);
								records[R].setReadName(readName);
								records[R].setMappingQuality(60);
								records[R].setProperPairFlag(true);
								records[R].setMateReferenceName(ssr.getSequenceName());
								records[R].setMateUnmappedFlag(false);
								records[R].setReadUnmappedFlag(false);
								records[R].setAttribute("RG", sample);
								records[R].setReadNegativeStrandFlag(R==1);
								
							
							 
								StringBuilder sequence = new StringBuilder( this.readSize);
								int readLen=0;
								int refPos1=(R==0?x:x+this.fragmentSize-this.readSize);
								records[R].setAlignmentStart(refPos1);
								
								List<CigarElement> cigarElements = new ArrayList<>( this.readSize);
								int NM=0;
								while(readLen< this.readSize)
									{
									String base=null;
									VariantContext overlap=null;
									for(VariantContext vc:variantBuffer)
										{
										int d= vc.getStart() - ( refPos1 + readLen) ;
										if( d==0)
											{
											overlap=vc;
											break;
											}
										if( d> 0) break;
										}
									
									if(overlap!=null)
										{
										Genotype genotype = overlap.getGenotype(sample);
										if(genotype.isCalled() && !genotype.isMixed())
											{
											List<Allele> alleles = genotype.getAlleles();
											if(alleles.size()!=2) throw new RuntimeException("Not a diploid organism.");
											Allele allele=null;
											if(genotype.isPhased())
												{
												allele =alleles.get(g_strand);
												}
											else
												{
												allele = alleles.get(Math.random()< 0.5?0:1);
												}
											if(allele.isSymbolic()) throw new  RuntimeException("Cannot handle symbolic alleles");
											if(allele.isReference())
												{
												cigarElements.add(new CigarElement(allele.getBaseString().length(), CigarOperator.EQ));
												}
											else if(allele.getBaseString().length() < overlap.getReference().getBaseString().length())
												{
												cigarElements.add(new CigarElement(allele.getBaseString().length(), CigarOperator.D));
												NM++;
												}
											else if(allele.getBaseString().length() > overlap.getReference().getBaseString().length())
												{
												cigarElements.add(new CigarElement(allele.getBaseString().length(), CigarOperator.I));
												NM++;
												}
											else //same size
												{
												cigarElements.add(new CigarElement(allele.getBaseString().length(), CigarOperator.X));
												NM++;
												}
											base=allele.getBaseString().toLowerCase();
											}
										}
									
									if(base==null)
										{
										base=String.valueOf(genomicSequence.charAt(refPos1-1+readLen));
										cigarElements.add(new CigarElement(1, CigarOperator.EQ));
										}
									sequence.append(base);
									readLen+=base.length();
									}//end loop over read-leng
								
								while(sequence.length()> this.readSize)
									sequence.deleteCharAt(sequence.length()-1);
								records[R].setReadString(sequence.toString());
								
								
								StringBuilder qual=new StringBuilder(sequence.length());
								for(int i=0;i< sequence.length();++i) qual.append("I");
								records[R].setBaseQualityString(qual.toString());
								
								//cigar
								int c=0;
								while(c+1< cigarElements.size())
									{
									CigarElement c0=cigarElements.get(c);
									CigarElement c1=cigarElements.get(c+1);
									if(c0.getOperator().equals(c1.getOperator()))
										{
										cigarElements.set(c, new CigarElement(c0.getLength()+c1.getLength(), c0.getOperator()));
										cigarElements.remove(c+1);
										}
									else
										{
										++c;
										}
									}
								records[R].setCigar(new Cigar(cigarElements));
								records[R].setAttribute("NM",NM);
								}//end loop over R1/R2
							
							if(Math.random()<0.5)
								{	 
								records[0].setFirstOfPairFlag(true);records[0].setSecondOfPairFlag(false);
								records[1].setFirstOfPairFlag(false);records[1].setSecondOfPairFlag(true);
								}
							else
								{
								records[0].setFirstOfPairFlag(false);records[0].setSecondOfPairFlag(true);
								records[1].setFirstOfPairFlag(true);records[1].setSecondOfPairFlag(false);
								}
							
							records[0].setMateAlignmentStart( records[1].getAlignmentStart() );
							records[1].setMateAlignmentStart( records[0].getAlignmentStart() );
							
							records[0].setInferredInsertSize(records[1].getAlignmentStart() - records[0].getAlignmentStart());
							records[1].setInferredInsertSize(records[0].getAlignmentStart() - records[1].getAlignmentStart());
		
							
							samFileWriter.addAlignment(records[0]);
							samFileWriter.addAlignment(records[1]);
							}
						}
					}
				
				++x;
				}
			
			}
		
		samFileWriter.close();
		}
	
	
	
	
	@Override
	public int doWork(final List<String> args)
		{
		if(this.faidx==null)
			{
			LOG.error("No REF defined");
			return -1;
			}
		if(readSize<=0) {
			LOG.error("bad read size");
			return -1;
			}
		if(fragmentSize<=readSize*2) {
			LOG.error("bad fragment size");
			return -1;
			}
		VCFIterator iter=null;
		try
			{
			this.indexedFastaSequenceFile= new IndexedFastaSequenceFile(this.faidx);
			iter  = super.openVCFIterator(super.oneFileOrNull(args));
			run(iter);
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			}
		}

	
	
	public static void main(String[] args)
		{
		new VcfToBam().instanceMainWithExit(args);
		}

	}
