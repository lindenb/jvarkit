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
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.biostar;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.io.PrintWriter;
import htsjdk.samtools.util.CloserUtil;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.command.Command;

public class Biostar94573 extends AbstractBiostar94573
	{
	private static final char CLIPPING=' ';
	private static final char DELETION='-';
	private static final char MATCH='*';
	private enum Format{None,Clustal,Fasta};
	
	
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(Biostar94573.class);

	@Override
	public Command createCommand() {
		return new MyCommand();
		}

	static private class MyCommand extends AbstractBiostar94573.AbstractBiostar94573Command
		{    
	
		private int align_length=0;
		private Map<String, AlignSequence> sample2sequence=new HashMap<String, AlignSequence>();
		private AbstractSequence consensus=null;
		
		
		private abstract class AbstractSequence
			{
			abstract char at(int index);
			@Override
			public String toString()
				{
				StringBuilder b=new StringBuilder(align_length);
				for(int i=0;i< align_length;++i) b.append(at(i));
				return b.toString();
				}
			}
		
		
		
		private abstract class Sequence extends AbstractSequence
			{
			StringBuilder seq=new StringBuilder();
			char at(int index)
				{
				return(index< 0 || index >=seq.length()?CLIPPING:Character.toUpperCase(seq.charAt(index)));
				}
			}
		
		private class AlignSequence extends Sequence
			{
			String name;
			}
		
		private interface Consensus
			{
			
			}
		
		private class ClustalConsensus extends Sequence implements Consensus
			{
			}
		
	
		private class FastaConsensus extends AbstractSequence implements Consensus
			{
			@Override
			char at(int index)
				{
				Counter<Character> count=new Counter<Character>();
				for(AlignSequence a:sample2sequence.values()) count.incr(a.at(index));
				if(count.getCountCategories()<=1) return MATCH;
				return ' ';
				}
			}
			
		private class NamedConsensus extends AbstractSequence implements Consensus
			{
			private AlignSequence namedConsensus;
			NamedConsensus(AlignSequence namedConsensus)
				{
				this.namedConsensus=namedConsensus;
				}
			@Override
			char at(int index)
				{
				for(AlignSequence a:sample2sequence.values())
					{
					if(a==this.namedConsensus) continue;
					if( a.at(index)!=this.namedConsensus.at(index)) return ' ';
					}
				return MATCH;
				}
			}
		
		
		@Override
		protected Collection<Throwable> call(String inputName) throws Exception {
				VariantContextWriter w=null;
				LineIterator r=null;
				try
					{
					if(inputName==null)
						{
						LOG.info("Reading from stdin");
						r=IOUtils.openStreamForLineIterator(stdin());
						}
					else
						{
						LOG.info("Reading from "+inputName);
						r=IOUtils.openURIForLineIterator(inputName);
						}
				
					Format format=Format.None;
		
					/** try to guess format */
					while(r.hasNext() && format==Format.None)
						{
						String line=r.peek();
						if( line.trim().isEmpty()) { r.next(); continue;}
						if(line.startsWith("CLUSTAL"))
							{
							format=Format.Clustal;
							r.next();//consume
							break;
							}
						else if(line.startsWith(">"))
							{
							format=Format.Fasta;
							break;
							}
						else
							{
							return wrapException("MSA format not recognized in "+line);
							}
						}
					LOG.info("format : "+format);
					/** parse lines as FASTA */
					if(Format.Fasta.equals(format))
						{
						this.consensus=new FastaConsensus();
						AlignSequence curr=null;
						while(r.hasNext())
							{
							String line=r.next();
							if(line.startsWith(">"))
								{
								curr=new AlignSequence();
								curr.name=line.substring(1).trim();
								if(sample2sequence.containsKey(curr.name))
									{
									return wrapException("Sequence ID "+curr.name +" defined twice");
									}
								sample2sequence.put(curr.name, curr);
								}
							else if(curr!=null)
								{
								curr.seq.append(line.trim());
								this.align_length=Math.max(this.align_length, curr.seq.length());
								}
							}
						}
					/** parse lines as CLUSTAL */
					else if(Format.Clustal.equals(format))
						{
						ClustalConsensus clustalconsensus=new ClustalConsensus();
						this.consensus=clustalconsensus;
						AlignSequence curr=null;
						int columnStart=-1;
						while(r.hasNext())
							{
							String line=r.next();
							
							if( line.trim().isEmpty() || line.startsWith("CLUSTAL W"))
								{
								columnStart=-1;
								continue;
								}
							if(line.charAt(0)==' ')
								{
								if(columnStart==-1)
									{
									return wrapException("illegal consensus line for "+line);
									}	
								/* if consensus doesn't exist in the first rows */
								while(clustalconsensus.seq.length() < (this.align_length-(line.length()-columnStart) ))
									{
									clustalconsensus.seq.append(" ");
									}
								clustalconsensus.seq.append(line.substring(columnStart));
								}
							else
								{
								 if(columnStart==-1)
									 {
									columnStart=line.indexOf(' ');
									if(columnStart==-1)
										{
										return wrapException("no whithespace in "+line);
										}
									while(columnStart< line.length() && line.charAt(columnStart)==' ')
										{
										columnStart++;
										}
									}
								String seqname=line.substring(0, columnStart).trim();
								curr=this.sample2sequence.get(seqname);
								if(curr==null)
									{
									curr=new AlignSequence();
									curr.name=seqname;
									this.sample2sequence.put(curr.name, curr);
									}
								int columnEnd=line.length();
								//remove blanks and digit at the end
								while(columnEnd-1>columnStart && (line.charAt(columnEnd-1)==' ' || Character.isDigit(line.charAt(columnEnd-1))))
										{
										columnEnd--;
										}
								curr.seq.append(line.substring(columnStart,columnEnd));
								this.align_length=Math.max(align_length, curr.seq.length());
								}
							}
						}
					else
						{
						return wrapException("Undefined input format");
						}
					CloserUtil.close(r);
					
					/* sequence consensus was set*/
					if(super.consensusRefName!=null)
						{	
						AlignSequence namedSequence=null;
						if((namedSequence=sample2sequence.get(super.consensusRefName))==null)
							{
							return wrapException("Cannot find consensus sequence \""+consensusRefName+"\" in list of sequences: "+this.sample2sequence.keySet().toString());
							}
						this.consensus = new NamedConsensus(namedSequence);
						}
					
					/** we're done, print VCF */
					
					/** first, print header */
					Set<VCFHeaderLine> vcfHeaderLines=new HashSet<VCFHeaderLine>();
		
					vcfHeaderLines.add(new VCFInfoHeaderLine(VCFConstants.DEPTH_KEY, 1, VCFHeaderLineType.Integer, "Approximate read depth."));
					vcfHeaderLines.add(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1, VCFHeaderLineType.String, "Genotype"));
					vcfHeaderLines.add(new VCFFormatHeaderLine(VCFConstants.DEPTH_KEY, 1, VCFHeaderLineType.Integer, "Approximate read depth"));
					vcfHeaderLines.add(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
					vcfHeaderLines.add(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
					Map<String,String> mapping=new HashMap<String,String>();
					mapping.put("ID", REF);
					mapping.put("length",String.valueOf(this.align_length));
					vcfHeaderLines.add(new VCFContigHeaderLine(mapping,1));
		
					Set<String> samples=new TreeSet<String>(this.sample2sequence.keySet());
					VCFHeader vcfHeader=new VCFHeader(vcfHeaderLines,samples);
					
					w = super.openVariantContextWriter(); 
					w.writeHeader(vcfHeader);
					
					/** loop over data, print header */
					int pos1=0;
					while(pos1< align_length)
						{
						boolean is_variation;//is it a real variation or print all sites
						if(consensus.at(pos1)==MATCH)
							{
							if(this.printAllSites)
								{
								is_variation=false;
								}
							else
								{
								++pos1;
								continue;
								}
							}
						else
							{
							is_variation=true;
							}
						int pos2=pos1+1;
						
						// don't extend if no variation and printAllSites
						while(is_variation && pos2<align_length && consensus.at(pos2)!=MATCH)
							{
							++pos2;
							}
							
						boolean is_subsitution=(pos1+1==pos2);
						if(is_subsitution && pos1!=0 && is_variation)//need pos1>0 because ALT contains prev base.
							{
							for(Sequence seq: this.sample2sequence.values())
								{
								if(seq.at(pos1)==DELETION)
									{
									is_subsitution=false;
									break;
									}
		 						}
							}
						
						Set<Allele> alleles=new HashSet<Allele>();
						
						VariantContextBuilder vcb=new VariantContextBuilder();
						List<Genotype> genotypes=new ArrayList<Genotype>(samples.size());
						
						/* longest variant */
						String longest=null;
						Counter<String> countAlleles=new Counter<String>();
						Map<String,String> sample2genotype=new HashMap<String,String>(samples.size());
						
						String namedConsensusRefAllele="N";
						
						/* loop over the sequences of each seample */
						for(String sample:samples)
							{
							Sequence seq=this.sample2sequence.get(sample);
							String al=null;
							if(is_subsitution)
								{
								if(seq.at(pos1)==CLIPPING) continue;
								al=String.valueOf(seq.at(pos1));
								}
							else
								{
								StringBuilder sb=new StringBuilder(pos2-pos1);
								for(int i=Math.max(0,pos1-1);//yes -1
										i<pos2;
										++i)
									{
									if(seq.at(i)==CLIPPING) continue;
									if(seq.at(i)==DELETION) continue;
									sb.append(seq.at(i));
									}
								if(sb.length()==0) continue;
								al=sb.toString();
								}
							/* did we find the longest allele ?*/
							if(longest==null || longest.length()< al.length())
								{
								countAlleles=new Counter<String>();//reset count of most frequent, we'll use the longest indel or subst
								longest=al;
								}
							countAlleles.incr(al);
							sample2genotype.put(sample,al);
							
							/* if consensus sequence name was defined , record this allele for future use */
							if(consensusRefName!=null && sample.equals(consensusRefName))
								{
								namedConsensusRefAllele = al;
								}
							
							}
						
						if(countAlleles.isEmpty())
							{
							if(printAllSites==false)//printAllSites=false
								{
								continue;
								}
							/* no a real variation, just add a dummy 'N' */
							countAlleles.incr("N");
							}
						String refAllStr;
						if( consensusRefName == null)
							{
							refAllStr = countAlleles.getMostFrequent();
							}
						else
							{
							refAllStr = namedConsensusRefAllele;
							}
						final Allele refAllele=Allele.create(refAllStr.replaceAll("[^ATGCatgc]","N"), true);
						alleles.add(refAllele);
						
						/* loop over samples, , build each genotype */
						for(String sample:sample2genotype.keySet())
							{
							
							Allele al=null;
							
							if(!sample2genotype.containsKey(sample))
								{
								//nothing
								}
							else if( sample2genotype.get(sample).equals(refAllStr))
								{
								al=refAllele;
								}
							else
								{
								al=Allele.create(sample2genotype.get(sample).replaceAll("[^ATGCatgc]","N"), false);
								alleles.add(al);
								}
							
							if(al!=null)
								{
								final GenotypeBuilder gb=new GenotypeBuilder(sample);
								final List<Allele> sampleAlleles=new ArrayList<Allele>(2);
								sampleAlleles.add(al);
								if(!haploid) sampleAlleles.add(al);
								gb.alleles(sampleAlleles);
								gb.DP(1);
								genotypes.add(gb.make());
								}
							else
								{
								genotypes.add(GenotypeBuilder.createMissing(sample, haploid?1:2));
								}
							
							
							}
						final int start=pos1+(is_subsitution?1:0);//got to 1-based ref if subst, for indel with use pos(base)-1
						vcb.start(start);
						vcb.stop(start+(refAllStr.length()-1));
						vcb.chr(REF);
						HashMap<String, Object> atts=new HashMap<String,Object>();
						atts.put(VCFConstants.DEPTH_KEY, genotypes.size());
						vcb.attributes(atts);
						vcb.alleles(alleles);
						vcb.genotypes(genotypes);
						w.add(vcb.make());
						pos1=pos2;
						}
					w.close();
					if(super.outFasta!=null)
						{
						PrintWriter fasta=new PrintWriter(outFasta);
						for(String sample:samples)
							{
							fasta.println(">"+sample);
							Sequence seq=this.sample2sequence.get(sample);
							for(int i=0;i< align_length;++i)
								{
								fasta.print(seq.at(i));
								}
							fasta.println();
							}	
						fasta.println(">CONSENSUS");
						for(int i=0;i< align_length;++i)
								{
								fasta.print(consensus.at(i));
								}
						fasta.println();
						fasta.flush();
						fasta.close();
						}
					
					LOG.info("Done");
					return Collections.emptyList();
					}
				catch(Exception err)
					{
					return wrapException(err);
					}
				finally
					{
					CloserUtil.close(r);
					CloserUtil.close(w);
					}
				}
			}
		
	public static void main(String[] args) {
		new Biostar94573().instanceMainWithExit(args);
	}
}
