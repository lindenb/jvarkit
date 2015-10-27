/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeSet;

import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class VcfToBam extends AbstractVcfToBam
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(VcfToBam.class);

	 @Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	 private static class MyCommand extends AbstractVcfToBam.AbstractVcfToBamCommand
	 	{
		private IndexedFastaSequenceFile indexedFastaSequenceFile = null;
		private int fragmentSize=600;
		private int readSize=100;
		
		
		private void run(VcfIterator vcfIterator) throws IOException
			{
			 long id_generator=0L;
			
			
			VCFHeader header = vcfIterator.getHeader();
			SAMSequenceDictionary dict = header.getSequenceDictionary();
			if(dict==null) throw new IOException("Sequence Dictionary missing in VCF");
			if(!SequenceUtil.areSequenceDictionariesEqual(dict, indexedFastaSequenceFile.getSequenceDictionary()))
				{
				LOG.warn("REF/VCF: not the same Sequence dictonary "+dict+" "+ indexedFastaSequenceFile.getSequenceDictionary());
				}
			SAMFileHeader samHeader=new SAMFileHeader();
			samHeader.setSequenceDictionary(dict);
			for(String sample:new TreeSet<>(header.getSampleNamesInOrder()))
				{
				SAMReadGroupRecord rg= new SAMReadGroupRecord(sample);
				rg.setSample(sample);
				rg.setLibrary(sample);
				rg.setDescription(sample);
				rg.setLibrary("illumina");
				samHeader.addReadGroup(rg);
				
				}
			samHeader.addComment("Generated with "+getProgramCommandLine());
			samHeader.setSortOrder(SortOrder.unsorted);
			SAMFileWriter samFileWriter= super.openSAMFileWriter(samHeader, true);
			
			/* looping over sequences */
			for(SAMSequenceRecord ssr: dict.getSequences())
				{
				LinkedList<VariantContext> variantBuffer = new LinkedList<>();
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
						
						VariantContext ctx = vcfIterator.peek();
						if(ctx==null) break;
						if(ctx.isIndel() || !ctx.isVariant())
							{
							LOG.warn("Cannot handle "+ctx.getReference()+"/"+ctx.getAlternateAlleles());
							vcfIterator.next();//consumme
							continue;
							}
						
						SAMSequenceRecord variantContig = dict.getSequence( ctx.getContig());
						if(variantContig == null) throw new IOException("Unknown contig "+ctx.getContig());
						if(variantContig.getSequenceIndex()< ssr.getSequenceIndex())
							{
							throw new IOException("Variants are not ordered on sequence dictionary ("+ctx+")");
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
						LOG.warn("no variant ?");
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
													}
												else if(allele.getBaseString().length() > overlap.getReference().getBaseString().length())
													{
													cigarElements.add(new CigarElement(allele.getBaseString().length(), CigarOperator.I));
													}
												else //same size
													{
													cigarElements.add(new CigarElement(allele.getBaseString().length(), CigarOperator.X));
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
		public Collection<Throwable> initializeKnime()
			{
			if(this.reference==null)
				{
				return wrapException("No REF defined");
				}
			try
				{
				this.indexedFastaSequenceFile = new IndexedFastaSequenceFile(this.reference);
				}
			catch(Exception err)
				{
				return wrapException(err);
				}
			return super.initializeKnime();
			}
		@Override
		public void disposeKnime() {
			CloserUtil.close(this.indexedFastaSequenceFile);
			this.indexedFastaSequenceFile=null;
			super.disposeKnime();
			}
		@Override
		protected Collection<Throwable> call(String inputName) throws Exception {
			VcfIterator iter=null;
			try
				{
				iter =  openVcfIterator(inputName);
				run(iter);
				return RETURN_OK;
				}
			catch(Exception err)
				{
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(iter);
				}
			}
	 	}
	
	
	public static void main(String[] args)
		{
		new VcfToBam().instanceMainWithExit(args);
		}

	}
