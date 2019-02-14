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
package com.github.lindenb.jvarkit.tools.msa2vcf;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.io.File;
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

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**

BEGIN_DOC


Deprecated: use https://github.com/sanger-pathogens/snp_sites , though some people told me they still use it for misc reasons.


## Motivation

Getting a VCF file from a CLUSTAW alignment. See also http://www.biostars.org/p/94573/

input is a clustalw file like: https://github.com/biopython/biopython/blob/master/Tests/Clustalw/opuntia.aln


## Cited-in


  * 'Differential distribution of Neandertal genomic signatures in human mitochondrial haplogroups'. 2017. Renata C Ferreira, Camila R Rodrigues, James R Broach, View ORCID ProfileMarcelo RS Briones. doi: [https://doi.org/10.1101/190363]([https://doi.org/10.1101/190363)
  * 'Pleiotropic effects of regulatory variation in tan result in correlation of two pigmentation traits in Drosophila melanogaster'. 2018. Molecular Ecology. Lukas Endler, Jean‐Michel Gibert, Viola Nolte, Christian Schlotterer. doi: 10.1111/mec.14781
  * 'Two key events associated with a transposable element burst occurred during rice domestication'  https://doi.org/10.1101/405290  
  * 'Tracking the origin of two genetic components associated with transposable element bursts in domesticated rice'. Nature Communicationsvolume 10, Article number: 641 (2019) https://www.nature.com/articles/s41467-019-08451-3

## Example

```bash
$ curl https://raw.github.com/biopython/biopython/master/Tests/Clustalw/opuntia.aln

CLUSTAL W (1.81) multiple sequence alignment


gi|6273285|gb|AF191659.1|AF191      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
gi|6273284|gb|AF191658.1|AF191      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
gi|6273287|gb|AF191661.1|AF191      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
gi|6273286|gb|AF191660.1|AF191      TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
gi|6273290|gb|AF191664.1|AF191      TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
gi|6273289|gb|AF191663.1|AF191      TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
gi|6273291|gb|AF191665.1|AF191      TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
                                    ******* **** *************************************

gi|6273285|gb|AF191659.1|AF191      TATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATA
gi|6273284|gb|AF191658.1|AF191      TATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATA
gi|6273287|gb|AF191661.1|AF191      TATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATA
gi|6273286|gb|AF191660.1|AF191      TATATA----------ATATATTTATAATTTCCTTATATATCCAAATATA
gi|6273290|gb|AF191664.1|AF191      TATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
gi|6273289|gb|AF191663.1|AF191      TATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
gi|6273291|gb|AF191665.1|AF191      TATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATA
                                    ******          ********  **** ********* *********

gi|6273285|gb|AF191659.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGT
gi|6273284|gb|AF191658.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGT
gi|6273287|gb|AF191661.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGT
gi|6273286|gb|AF191660.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGT
gi|6273290|gb|AF191664.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGT
gi|6273289|gb|AF191663.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTAT
gi|6273291|gb|AF191665.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGT
                                    ************************************ *********** *

gi|6273285|gb|AF191659.1|AF191      ACCAGA
gi|6273284|gb|AF191658.1|AF191      ACCAGA
gi|6273287|gb|AF191661.1|AF191      ACCAGA
gi|6273286|gb|AF191660.1|AF191      ACCAGA
gi|6273290|gb|AF191664.1|AF191      ACCAGA
gi|6273289|gb|AF191663.1|AF191      ACCAGA
gi|6273291|gb|AF191665.1|AF191      ACCAGA
                                    ******
```
generate the VCF

```
$ curl https://raw.github.com/biopython/biopython/master/Tests/Clustalw/opuntia.aln" |\
  java -jar dist/msa2vcf.jar

##fileformat=VCFv4.1
##Biostar94573CmdLine=
##Biostar94573Version=ca765415946f3ed0827af0773128178bc6aa2f62
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth.">
##contig=<ID=chrUn,length=156>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	gi|6273284|gb|AF191658.1|AF191	gi|6273285|gb|AF191659.1|AF191	gi|6273286|gb|AF191660.1|AF191	gi|6273287|gb|AF191661.1|AF191	gi|6273289|gb|AF191663.1|AF191	gi|6273290|gb|AF191664.1|AF191	gi|6273291|gb|AF191665.1|AF191
chrUn	8	.	T	A	.	.	DP=7	GT:DP	0:1	0:1	1:1	0:1	0:1	0:1	0:1
chrUn	13	.	A	G	.	.	DP=7	GT:DP	0:1	0:1	0:1	0:1	1:1	1:1	1:1
chrUn	56	.	ATATATATATA	ATA,A,ATATA	.	.	DP=7	GT:DP	1:1	2:1	2:1	2:1	3:1	3:1	0:1
chrUn	74	.	TCA	TAT	.	.	DP=7	GT:DP	0:1	0:1	1:1	0:1	0:1	0:1	0:1
chrUn	81	.	T	C	.	.	DP=7	GT:DP	0:1	0:1	0:1	0:1	1:1	1:1	1:1
chrUn	91	.	T	C	.	.	DP=7	GT:DP	1:1	1:1	0:1	0:1	0:1	0:1	0:1
chrUn	137	.	T	C	.	.	DP=7	GT:DP	0:1	1:1	0:1	0:1	0:1	0:1	0:1
chrUn	149	.	G	A	.	.	DP=7	GT:DP	0:1	0:1	0:1	0:1	1:1	0:1	0:1
```

END_DOC
*/


@Program(name="msa2vcf",
	description="Getting a VCF file from a CLUSTAW or a FASTA alignment. ",
	deprecatedMsg="use https://github.com/sanger-pathogens/snp_sites",
	biostars=94573,
	keywords={"vcf","snp","msa","alignment"}
	)
public class MsaToVcf extends Launcher
	{
	private static final Logger LOG = Logger.build(MsaToVcf.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-R","--REF"},description="reference name used for the CHROM column. Optional")
	private String REF = "chrUn";

	@Parameter(names={"-c","--consensus"},description="ruse this sequence as CONSENSUS")
	private String consensusRefName = null;

	@Parameter(names={"-f","--fasta"},description="save computed fasta sequence in this file.")
	private File outFasta = null;

	@Parameter(names={"-m","--haploid"},description="haploid output")
	private boolean haploid = false;

	@Parameter(names={"-a","--allsites"},description="print all sites")
	private boolean printAllSites = false;

	private static final char CLIPPING=' ';
	private static final char DELETION='-';
	private static final char MATCH='*';
	private int align_length=0;
	private Map<String, AlignSequence> sample2sequence=new HashMap<String, AlignSequence>();
	private AbstractSequence consensus=null;
	private enum Format{None,Clustal,Fasta};
	
	
	private abstract class AbstractSequence
		{
		abstract char at(int index);
		@Override
		public String toString()
			{
			final StringBuilder b=new StringBuilder(align_length);
			for(int i=0;i< align_length;++i) b.append(at(i));
			return b.toString();
			}
		}
	
	
	
	private abstract class Sequence extends AbstractSequence
		{
		StringBuilder seq=new StringBuilder();
		@Override char at(int index)
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
			final Counter<Character> count=new Counter<Character>();
			for(final AlignSequence a:sample2sequence.values()) count.incr(a.at(index));
			if(count.getCountCategories()<=1) return MATCH;
			return ' ';
			}
		}
		
	private class NamedConsensus extends AbstractSequence implements Consensus
		{
		private final AlignSequence namedConsensus;
		NamedConsensus(final AlignSequence namedConsensus)
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
	public int doWork(final List<String> args) {
		VariantContextWriter w=null;
		LineIterator r=null;
		try
			{
			final String inputName= oneFileOrNull(args);
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
				final String line=r.peek();
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
					LOG.error("MSA format not recognized in "+line);
					return -1;
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
							LOG.error("Sequence ID "+curr.name +" defined twice");
							return -1;
							}
						sample2sequence.put(curr.name, curr);
						}
					else if(curr!=null)
						{
						curr.seq.append(line.trim());
						this.align_length=Math.max(this.align_length, curr.seq.length());
						}
					}
				/*
				//remove heading & trailing '-'
				for(final String sample:this.sample2sequence.keySet())
					{
					final AlignSequence seq = this.sample2sequence.get(sample);
					int i=0;
					while(i<this.align_length && seq.at(i)==DELETION)
						{
						seq.seq.setCharAt(i, CLIPPING);
						++i;
						}
					i= this.align_length-1;
					while(i>=0 && seq.at(i)==DELETION)
						{
						seq.seq.setCharAt(i, CLIPPING);
						--i;
						}
					}*/
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
							LOG.error("illegal consensus line for "+line);
							return -1;
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
								LOG.error("no whithespace in "+line);
								return -1;
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
				LOG.error("Undefined input format");
				return -1;
				}
			CloserUtil.close(r);
			
			/* sequence consensus was set*/
			if( consensusRefName!=null)
				{	
				AlignSequence namedSequence=null;
				if((namedSequence=sample2sequence.get(consensusRefName))==null)
					{
					LOG.error("Cannot find consensus sequence \""+consensusRefName+"\" in list of sequences: "+this.sample2sequence.keySet().toString());
					return -1;
					}
				this.consensus = new NamedConsensus(namedSequence);
				}
			/** we're done, print VCF */
			
			/** first, print header */
			Set<VCFHeaderLine> vcfHeaderLines=new HashSet<VCFHeaderLine>();

			vcfHeaderLines.add(new VCFInfoHeaderLine(VCFConstants.DEPTH_KEY, 1, VCFHeaderLineType.Integer, "Approximate read depth."));
			vcfHeaderLines.add(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1, VCFHeaderLineType.String, "Genotype"));
			vcfHeaderLines.add(new VCFFormatHeaderLine(VCFConstants.DEPTH_KEY, 1, VCFHeaderLineType.Integer, "Approximate read depth"));
			//super.addMetaData(vcfHeaderLines);
			Map<String,String> mapping=new HashMap<String,String>();
			mapping.put("ID", REF);
			mapping.put("length",String.valueOf(this.align_length));
			vcfHeaderLines.add(new VCFContigHeaderLine(mapping,1));

			Set<String> samples=new TreeSet<String>(this.sample2sequence.keySet());
			VCFHeader vcfHeader=new VCFHeader(vcfHeaderLines,samples);
			
			w= super.openVariantContextWriter(this.outputFile);
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
			if(outFasta!=null)
				{
				final PrintWriter fasta= super.openFileOrStdoutAsPrintWriter(outFasta);
				for(final String sample:samples)
					{
					fasta.println(">"+sample);
					final Sequence seq=this.sample2sequence.get(sample);
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
			
			
			
			return RETURN_OK;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);
			CloserUtil.close(w);
			}
		}
	
	public static void main(String[] args) {
		new MsaToVcf().instanceMainWithExit(args);
	}
}
