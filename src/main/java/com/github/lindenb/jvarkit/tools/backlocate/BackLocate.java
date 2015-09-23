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
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;



import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;
import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;

public class BackLocate
	extends AbstractBackLocate
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(BackLocate.class);

	
	@Override
	public Command createCommand()
		{
		return new MyCommand();
		}
	
	static private class MyCommand extends AbstractBackLocate.AbstractBackLocateCommand
		{
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
		PrintStream out,
		KnownGene gene,
		String geneName,
		char aa1,char aa2,
		int peptidePos1
		) throws IOException
		{
		
		final GeneticCode geneticCode=getGeneticCodeByChromosome(gene.getChromosome());
		RNASequence wildRNA=null;
		ProteinCharSequence wildProt=null;
		
	        		
	        		
		if(genomicSeq==null ||
	               !gene.getChromosome().equals(genomicSeq.getChrom()) 
	               )
        	{
        	LOG.info("fetch genome");
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
	    	 stderr().println("#no protein found for transcript:"+gene.getName());
	    	 return;
	    	 }
	    int peptideIndex0= peptidePos1-1;
        if(peptideIndex0 >=wildProt.length())
        	{
        	wrapException("#index out of range for :"+gene.getName()+" petide length="+wildProt.length());
        	}
    
        if(wildProt.charAt(peptideIndex0)!=aa1)
        	{
        	out.println("##Warning ref aminod acid for "+gene.getName() +"  ["+peptidePos1+"] is not the same ("+wildProt.charAt(peptideIndex0)+"/"+aa1+")");
        	}
        else
        	{
        	out.println("##"+gene.getName());
        	}
        int indexesInRNA[]=new int[]{
        	0+ peptideIndex0*3,
        	1+ peptideIndex0*3,
        	2+ peptideIndex0*3
        	};
        final String wildCodon=""
        		+ wildRNA.charAt(indexesInRNA[0])
        		+ wildRNA.charAt(indexesInRNA[1])
        		+ wildRNA.charAt(indexesInRNA[2])
        		;
        /* 2015 : adding possible mut codons */
        final Set<String> possibleAltCodons = new HashSet<>();
        final char bases[]=new char[]{'A','C','G','T'};
        for(int codon_pos=0;codon_pos<3;++codon_pos)
        	{
        	StringBuilder sb=new StringBuilder(wildCodon);
        	for(char mutBase:bases)
        		{
        		sb.setCharAt(codon_pos, mutBase);
        		if(geneticCode.translate(sb.charAt(0), sb.charAt(1), sb.charAt(2))==Character.toUpperCase(aa2))
        			{
        			possibleAltCodons.add(sb.toString());
        			}
        		}
        	}
        
        for(int indexInRna: indexesInRNA)
        	{
        	out.print(geneName);
        	out.print('\t');
        	out.print(aa1);
        	out.print('\t');
        	out.print(peptidePos1);
        	out.print('\t');
        	out.print(aa2);
        	out.print('\t');
        	out.print(gene.getName());
        	out.print('\t');
        	out.print(gene.getStrand()==Strand.NEGATIVE?"-":"+");
        	out.print('\t');
        	out.print(wildProt.charAt(peptideIndex0));
        	out.print('\t');
        	out.print(indexInRna);
        	out.print('\t');
        	out.print(wildCodon);
        	out.print('\t');
        	if(possibleAltCodons.isEmpty())
        		{
        		out.print('.');
        		}
        	else
        		{
        		boolean first=true;
        		for(String mutCodon:possibleAltCodons)
        			{
        			if(!first) out.print('|');
        			first=false;
        			out.print(mutCodon);
        			}
        		}
        	out.print('\t');
        	out.print(wildRNA.charAt(indexInRna));
        	out.print('\t');
        	out.print(gene.getChromosome());
        	out.print('\t');
        	out.print(wildRNA.genomicPositions.get(indexInRna));
        	out.print('\t');
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
        	out.print(exonName);
        	if(this.printSequences)
        		{
        		String s=wildRNA.toString();
        		out.print('\t');
            	out.print(s.substring(0,indexInRna)+"["+s.charAt(indexInRna)+"]"+(indexInRna+1<s.length()?s.substring(indexInRna+1):""));
            	s=wildProt.toString();
            	out.print('\t');
            	out.print(s.substring(0,peptideIndex0)+"["+aa1+"/"+aa2+"/"+wildProt.charAt(peptideIndex0)+"]"+(peptideIndex0+1<s.length()?s.substring(peptideIndex0+1):""));
        		}
        	out.println();
        	}
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
	
	private void run(PrintStream out,LineIterator in) throws IOException
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
				LOG.warn("No kgXref found for "+geneName);
				continue;
				}
			
			
			for(String kgId:kgIds)
				{
				KnownGene kg=this.knwonGenes.get(kgId);
				if(kg==null) continue;
				backLocate(out,kg, geneName, aa1, aa2, position1);
				}
			}
		}
	
	
	private void loadKnownGenesFromUri(String kgURI) throws IOException
		{
		if(this.indexedFastaSequenceFile.getSequenceDictionary()==null)
			{
			throw new IOException("Cannot get sequence dictionary for REF : "+getMessageBundle("picard.dictionary.needed"));
			}
		
		LOG.info("loading genes");
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
					LOG.warn("The reference doesn't contain chromosome "+g.getContig());
					unknown.add(g.getContig());
					}
				continue;
				}
			
			this.knwonGenes.put(g.getName(),g);
			}
		in.close();
		LOG.info("genes:"+this.knwonGenes.size());
		}
	
	private void loadkgXRefFromUri(String kgURI) throws IOException
		{
		
		LOG.info("loading "+kgURI);
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
		LOG.info("kgxref:"+geneSymbol2kg.size());
		}

	
	@Override
	public Collection<Throwable> initializeKnime()
		{
		try
			{
			if(getReferenceFile()==null)
				{
				return wrapException(getMessageBundle("reference.undefined"));
				}
			this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(getReferenceFile());
			
			if(knownGeneURI==null)
				{
				return wrapException("Undefined knwonGeneURI");
				}
			
			if(kgXRef==null)
				{
				return wrapException("Undefined kgXref");
				}
			this.loadKnownGenesFromUri(knownGeneURI);
			this.loadkgXRefFromUri(kgXRef);
			
			return super.initializeKnime();
			}
		catch (Exception e)
			{
			return wrapException(e);
			}
		
		}
	
	@Override
	public Collection<Throwable> call() throws Exception
		{
		final List<String> args = this.getInputFiles();
		PrintStream out=null;
		try {			
			
			
			out = this.openFileOrStdoutAsPrintStream();
			
			out.print("#User.Gene");
        	out.print('\t');
        	out.print("AA1");
        	out.print('\t');
        	out.print("petide.pos.1");
        	out.print('\t');
        	out.print("AA2");
        	out.print('\t');
        	out.print("knownGene.name");
        	out.print('\t');
        	out.print("knownGene.strand");
        	out.print('\t');
        	out.print("knownGene.AA");
        	out.print('\t');
        	out.print("index0.in.rna");
        	out.print('\t');
        	out.print("wild.codon");
        	out.print('\t');
        	out.print("potential.var.codons");
        	out.print('\t');
        	out.print("base.in.rna");
        	out.print('\t');
        	out.print("chromosome");
        	out.print('\t');
        	out.print("index0.in.genomic");
        	out.print('\t');
        	out.print("exon");
        	if(this.printSequences)
        		{
        		out.print('\t');
            	out.print("mRNA");
            	out.print('\t');
            	out.print("protein");
        		}
        	out.println();
			if(args.isEmpty())
				{
				LOG.info("reading from stdin");
				LineIterator in=IOUtils.openStdinForLineIterator();
				this.run(out,in);
				CloserUtil.close(in);
				}
			else
				{
				for(String filename:args)
					{
					LOG.info("reading from "+filename);
					LineIterator in=IOUtils.openURIForLineIterator(filename);
					this.run(out,in);
					CloserUtil.close(in);
					}
				}
			return Collections.emptyList();
			}
		catch (Exception e) {
			return wrapException(e);
			}
		finally
			{
			CloserUtil.close(out);
			}	
		}
		
		@Override
		public void disposeKnime()
			{
			CloserUtil.close(this.indexedFastaSequenceFile);
			this.indexedFastaSequenceFile=null;
			super.disposeKnime();
			}
		
		}
	public static void main(String[] args)
		{
		new BackLocate().instanceMainWithExit(args);
		}
	}
