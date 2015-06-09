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
* Nov 2014: removed dependencies to SQL

*/
package com.github.lindenb.jvarkit.tools.backlocate;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.readers.LineIterator;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;



import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;

public class BackLocate
	extends AbstractCommandLineProgram
	{
	private static final String DEFAULT_KGXREF= "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/kgXref.txt.gz";
	private static final String DEFAULT_KNOWNGENE= "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz";
	
	private boolean printSequences=false;
	private GenomicSequence genomicSeq=null;
	private Map<String,Set<String>> geneSymbol2kg=new HashMap<>();
	private Map<String,KnownGene> knwonGenes=new HashMap<>();
	private IndexedFastaSequenceFile indexedFastaSequenceFile;
	
	/** get a genetic code from a chromosome name (either std or mitochondrial */
	private static GeneticCode getGeneticCodeByChromosome(String chr)
		{
		if(chr.equalsIgnoreCase("chrM") || chr.equalsIgnoreCase("MT")) return GeneticCode.getMitochondrial();
		return GeneticCode.getStandard();
		}




	
	static private class RNASequence extends AbstractCharSequence
		{
		List<Integer> genomicPositions=new ArrayList<Integer>();
		GenomicSequence genomic;
		char strand;
		RNASequence(GenomicSequence genomic,char strand)
			{
			this.genomic=genomic;
			this.strand=strand;
			}
		@Override
		public char charAt(int i)
			{
			char c=genomic.charAt(this.genomicPositions.get(i));
			return (strand=='+'?c:complement(c));
			}
		@Override
		public int length()
			{
			return genomicPositions.size();
			}
		}
	
	static private class ProteinCharSequence extends AbstractCharSequence
		{
		private RNASequence cDNA;
		private GeneticCode geneticCode;
		ProteinCharSequence(GeneticCode geneticCode,RNASequence cDNA)
			{
			this.geneticCode=geneticCode;
			this.cDNA=cDNA;
			}
		
		@Override
		public char charAt(int i)
			{
			return geneticCode.translate(
				cDNA.charAt(i*3+0),
				cDNA.charAt(i*3+1),
				cDNA.charAt(i*3+2));
			}	
		
		@Override
		public int length()
			{
			return this.cDNA.length()/3;
			}
	}


	

	

		

	private void backLocate(
		KnownGene gene,
		String geneName,
		char aa1,char aa2,
		int peptidePos1
		) throws IOException
		{
		
		GeneticCode geneticCode=getGeneticCodeByChromosome(gene.getChromosome());
		RNASequence wildRNA=null;
		ProteinCharSequence wildProt=null;
		
	        		
	        		
		if(genomicSeq==null ||
	               !gene.getChromosome().equals(genomicSeq.getChrom()) 
	               )
        	{
        	this.info("fetch genome");
        	this.genomicSeq= new GenomicSequence(this.indexedFastaSequenceFile, gene.getContig());
        	}
        	
	        		
	        		
	        		
	     if(gene.isPositiveStrand())
    		{    		
    		int exon_index=0;
    		while(exon_index< gene.getExonCount())
    			{
    			for(int i= gene.getExonStart(exon_index);
    					i< gene.getExonEnd(exon_index);
    					++i)
    				{
    				if(i< gene.getCdsStart()) continue;
    				if(i>=gene.getCdsEnd()) break;
					
					if(wildRNA==null)
						{
						wildRNA=new RNASequence(genomicSeq,'+');
						}

    				wildRNA.genomicPositions.add(i);
    				
    				
    				
    				if(wildRNA.length()%3==0 && wildRNA.length()>0 && wildProt==null)
        				{
        				wildProt=new ProteinCharSequence(geneticCode,wildRNA);
        				}
    				}
    			++exon_index;
    			}
    		
    		
    		
    		}
	   else // reverse orientation
    		{
    		int exon_index = gene.getExonCount()-1;
    		while(exon_index >=0)
    			{
    			for(int i= gene.getExonEnd(exon_index)-1;
    				    i>= gene.getExonStart(exon_index);
    				--i)
    				{
    				if(i>= gene.getCdsEnd()) continue;
    				if(i<  gene.getCdsStart()) break;
    				
    				if(wildRNA==null)
						{
						wildRNA=new RNASequence(genomicSeq,'-');
						}
    				
    				
    				
    				wildRNA.genomicPositions.add(i);
    				if( wildRNA.length()%3==0 &&
    					wildRNA.length()>0 &&
    					wildProt==null)
        				{
        				wildProt=new ProteinCharSequence(geneticCode,wildRNA);
        				}
    				
    				}
    			--exon_index;
    			}

    		}//end of if reverse
	        		
	     if(wildProt==null)
	    	 {
	    	 System.err.println("#no protein found for transcript:"+gene.getName());
	    	 return;
	    	 }
	    int peptideIndex0= peptidePos1-1;
        if(peptideIndex0 >=wildProt.length())
        	{
        	System.err.println("#index out of range for :"+gene.getName()+" petide length="+wildProt.length());
	    	return;
        	}
    
        if(wildProt.charAt(peptideIndex0)!=aa1)
        	{
        	System.out.println("##Warning ref aminod acid for "+gene.getName() +"  ["+peptidePos1+"] is not the same ("+wildProt.charAt(peptideIndex0)+"/"+aa1+")");
        	}
        else
        	{
        	System.out.println("##"+gene.getName());
        	}
        int indexesInRNA[]=new int[]{
        	0+ peptideIndex0*3,
        	1+ peptideIndex0*3,
        	2+ peptideIndex0*3
        	};
        String codon=""
        		+ wildRNA.charAt(indexesInRNA[0])
        		+ wildRNA.charAt(indexesInRNA[1])
        		+ wildRNA.charAt(indexesInRNA[2])
        		;
        		
        for(int indexInRna: indexesInRNA)
        	{
        	System.out.print(geneName);
        	System.out.print('\t');
        	System.out.print(aa1);
        	System.out.print('\t');
        	System.out.print(peptidePos1);
        	System.out.print('\t');
        	System.out.print(aa2);
        	System.out.print('\t');
        	System.out.print(gene.getName());
        	System.out.print('\t');
        	System.out.print(gene.getStrand()==Strand.NEGATIVE?"-":"+");
        	System.out.print('\t');
        	System.out.print(wildProt.charAt(peptideIndex0));
        	System.out.print('\t');
        	System.out.print(indexInRna);
        	System.out.print('\t');
        	System.out.print(codon);
        	System.out.print('\t');
        	System.out.print(wildRNA.charAt(indexInRna));
        	System.out.print('\t');
        	System.out.print(gene.getChromosome());
        	System.out.print('\t');
        	System.out.print(wildRNA.genomicPositions.get(indexInRna));
        	System.out.print('\t');
        	String exonName=null;
        	for(KnownGene.Exon exon : gene.getExons())
				{
				int genome=wildRNA.genomicPositions.get(indexInRna);
				if(exon.getStart()<=genome && genome< exon.getEnd())
					{
					exonName=exon.getName();
					break;
					}
				}
        	System.out.print(exonName);
        	if(this.printSequences)
        		{
        		String s=wildRNA.toString();
        		System.out.print('\t');
            	System.out.print(s.substring(0,indexInRna)+"["+s.charAt(indexInRna)+"]"+(indexInRna+1<s.length()?s.substring(indexInRna+1):""));
            	s=wildProt.toString();
            	System.out.print('\t');
            	System.out.print(s.substring(0,peptideIndex0)+"["+aa1+"/"+aa2+"/"+wildProt.charAt(peptideIndex0)+"]"+(peptideIndex0+1<s.length()?s.substring(peptideIndex0+1):""));
        		}
        	System.out.println();
        	}
		}

	private BackLocate() 
		{
		}
	
	private static char complement(char c)
		{
		switch(c)
			{
			case 'A':case 'a': return 'T';
			case 'T':case 't': return 'A';
			case 'G':case 'g': return 'C';
			case 'C':case 'c': return 'G';
			default:throw new IllegalArgumentException(""+c);
			}
		}
	
	private void run(LineIterator in) throws IOException
		{
		while(in.hasNext())
			{
			String line=in.next();
			if(line.startsWith("#") || line.trim().isEmpty()) continue;
			int n=line.indexOf('\t');
			if(n==0 || n==-1) throw new IOException("Bad line. No tab found in "+line);
			String geneName=line.substring(0,n).trim();
			if(geneName.isEmpty()) throw new IOException("Bad line. No gene in "+geneName);
			String mut=line.substring(n+1).trim();
			if(!mut.matches("[A-Za-z\\*][0-9]+[A-Za-z\\*]")) throw new IOException("Bad mutation  in "+line);
			char aa1= mut.substring(0,1).toUpperCase().charAt(0);
			char aa2= mut.substring(mut.length()-1).toUpperCase().charAt(0);
			int position1=Integer.parseInt(mut.substring(1,mut.length()-1));
			if(position1==0) throw new IOException("Bad position  in "+line);
			Set<String> kgIds= this.geneSymbol2kg.get(geneName.toUpperCase());
			if(kgIds==null || kgIds.isEmpty())
				{
				warning("No kgXref found for "+geneName);
				continue;
				}
			
			
			for(String kgId:kgIds)
				{
				KnownGene kg=this.knwonGenes.get(kgId);
				if(kg==null) continue;
				backLocate(kg, geneName, aa1, aa2, position1);
				}
			}
		}
	
	@Override
	protected String getOnlineDocUrl()
		{
		return "https://github.com/lindenb/jvarkit/wiki/BackLocate";
		}
	
	@Override
	public String getProgramDescription()
		{
		return "Mapping a mutation on a protein back to the genome.";
		}
	
	private void loadKnownGenesFromUri(String kgURI) throws IOException
		{
		if(this.indexedFastaSequenceFile.getSequenceDictionary()==null)
			{
			throw new IOException("Cannot get sequence dictionary for REF : "+getMessageBundle("picard.dictionary.needed"));
			}
		
		info("loading genes");
		Set<String> unknown=new HashSet<String>();
		BufferedReader in=IOUtils.openURIForBufferedReading(kgURI);
		String line;
		Pattern tab=Pattern.compile("[\t]");
		while((line=in.readLine())!=null)
			{
			if(line.isEmpty()) continue;
			String tokens[]=tab.split(line);
			KnownGene g=new KnownGene(tokens);
			Interval rgn=new Interval(g.getContig(), g.getTxStart()+1, g.getTxEnd());
			if(this.indexedFastaSequenceFile.getSequenceDictionary().getSequence(rgn.getContig())==null)
				{
				if(!unknown.contains(g.getContig()))
					{
					warning("The reference doesn't contain chromosome "+g.getContig());
					unknown.add(g.getContig());
					}
				continue;
				}
			
			this.knwonGenes.put(g.getName(),g);
			}
		in.close();
		info("genes:"+this.knwonGenes.size());
		}
	
	private void loadkgXRefFromUri(String kgURI) throws IOException
		{
		
		info("loading "+kgURI);
		BufferedReader in=IOUtils.openURIForBufferedReading(kgURI);
		String line;
		Pattern tab=Pattern.compile("[\t]");
		while((line=in.readLine())!=null)
			{
			if(line.isEmpty()) continue;
			String tokens[]=tab.split(line);
			String kgId=tokens[0];
			if(!this.knwonGenes.containsKey(kgId)) continue;
			String geneSymbol=tokens[4];
			Set<String> kglist= geneSymbol2kg.get(geneSymbol.toUpperCase());
			if(kglist==null)
				{
				kglist=new HashSet<String>();
				geneSymbol2kg.put(geneSymbol.toUpperCase(),kglist);
				}
			kglist.add(kgId);//kgID
			}
		in.close();
		info("kgxref:"+geneSymbol2kg.size());
		}

	@Override
	public void printOptions(PrintStream out)
		{
		System.out.println(" -R (fasta) "+getMessageBundle("reference.faidx"));
		System.out.println(" -k (knownGene) "+getMessageBundle("known.genes.uri")+" chromosomes must be named the same way than the REF sequence. Default "+DEFAULT_KGXREF);
		System.out.println(" -x (kgXRef uri)UCSC kgXref URI  Default:"+DEFAULT_KGXREF);
		System.out.println(" -p print mRNA & protein sequences");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		String knownGeneURI=null;
		String kgXref=null;
		
		try {			
			com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
			int c;
			while((c=opt.getopt(args,getGetOptDefault()+ "k:x:R:p"))!=-1)
				{
				switch(c)
					{
					case 'k': knownGeneURI=opt.getOptArg();break;
					case 'x': kgXref=opt.getOptArg();break;
					case 'R': this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(new File(opt.getOptArg()));break;
					case 'p': this.printSequences=true;break;
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
			if(this.indexedFastaSequenceFile==null)
				{
				error(getMessageBundle("reference.undefined"));
				return -1;
				}
			if(knownGeneURI==null)
				{
				warning("Undefined knwonGeneURI, using "+DEFAULT_KNOWNGENE);
				knownGeneURI=DEFAULT_KNOWNGENE;
				}
			
			if(kgXref==null)
				{
				warning("Undefined kgXref, using "+DEFAULT_KGXREF);
				kgXref=DEFAULT_KGXREF;
				}
			loadKnownGenesFromUri(knownGeneURI);
			loadkgXRefFromUri(kgXref);
			
			System.out.print("#User.Gene");
        	System.out.print('\t');
        	System.out.print("AA1");
        	System.out.print('\t');
        	System.out.print("petide.pos.1");
        	System.out.print('\t');
        	System.out.print("AA2");
        	System.out.print('\t');
        	System.out.print("knownGene.name");
        	System.out.print('\t');
        	System.out.print("knownGene.strand");
        	System.out.print('\t');
        	System.out.print("knownGene.AA");
        	System.out.print('\t');
        	System.out.print("index0.in.rna");
        	System.out.print('\t');
        	System.out.print("codon");
        	System.out.print('\t');
        	System.out.print("base.in.rna");
        	System.out.print('\t');
        	System.out.print("chromosome");
        	System.out.print('\t');
        	System.out.print("index0.in.genomic");
        	System.out.print('\t');
        	System.out.print("exon");
        	if(this.printSequences)
        		{
        		System.out.print('\t');
            	System.out.print("mRNA");
            	System.out.print('\t');
            	System.out.print("protein");
        		}
        	System.out.println();
			if(opt.getOptInd()==args.length)
				{
				info("reading from stdin");
				LineIterator in=IOUtils.openStdinForLineIterator();
				this.run(in);
				CloserUtil.close(in);
				}
			else
				{
				for(int optind=opt.getOptInd();optind<args.length;++optind)
					{
					String filename=args[optind++];
					info("reading from "+filename);
					LineIterator in=IOUtils.openURIForLineIterator(filename);
					this.run(in);
					CloserUtil.close(in);
					}
				}
			return 0;
			}
		catch (Exception e) {
			error(e);
			return -1;
			}
		finally
			{
			}	
		}
	public static void main(String[] args)
		{
		new BackLocate().instanceMainWithExit(args);
		}
	}
