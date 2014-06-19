package com.github.lindenb.jvarkit.tools.backlocate;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.readers.LineIterator;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.io.StringReader;
import java.net.URL;
import java.net.URLEncoder;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;

public class BackLocate
	extends AbstractCommandLineProgram
	{
	private boolean printSequences=false;
	private String genomeVersion="hg19";
	private Connection connection;
	private GenomicSequence genomicSeq=null;
	private DasSequenceProvider dasServer=new DasSequenceProvider();


	/** get a genetic code from a chromosome name (either std or mitochondrial */
	private static GeneticCode getGeneticCodeByChromosome(String chr)
		{
		if(chr.equalsIgnoreCase("chrM") || chr.equalsIgnoreCase("MT")) return GeneticCode.getMitochondrial();
		return GeneticCode.getStandard();
		}



	/**
	 * 
	 * A GenomicSequence
	 *
	 */
	static private class GenomicSequence
		extends AbstractCharSequence
		{
		private String chrom;
		private byte array[];
		private int chromStart0;
		
		public GenomicSequence(byte array[],String chrom,int chromStart0)
			{	
			this.chrom=chrom;
			this.array=array;
			this.chromStart0=chromStart0;
			}
		
		public String getChrom()
			{
			return chrom;
			}
		public int getChromStart()
			{
			return chromStart0;
			}
		public int getChromEnd()
			{
			return getChromStart()+array.length;
			}
		
		@Override
		public int length()
			{
			return getChromEnd();
			}
		
		@Override
		public char charAt(int index0)
			{
			if(index0 < getChromStart() || index0 >=getChromEnd())
				{
				throw new IndexOutOfBoundsException("index:"+index0);
				}
			return (char)(array[index0-chromStart0]);
			}
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


	/**
	 * Calls Ucsc DAS to fetch a DNA sequence using a SAX parser
	 */
	private class DasSequenceProvider
		extends DefaultHandler
		{
		private ByteArrayOutputStream baos=null;
		int reserve=100000;
		private SAXParser parser;
		public DasSequenceProvider()
			{
			SAXParserFactory f=SAXParserFactory.newInstance();
			f.setSchema(null);
			f.setNamespaceAware(false);
			f.setValidating(false);
			try
				{
				this.parser=f.newSAXParser();
				}
			catch(Exception err)
				{
				throw new RuntimeException(err);
				}
			}
		
		 @Override
		public void startDocument() throws SAXException
			{
			baos=null;
			}
		
		 public InputSource resolveEntity (String publicId, String systemId)
	         {
	         return new InputSource(new StringReader(""));
	         }
		
		@Override
	    public void startElement(String uri, String localName, String name,
	            Attributes attributes) throws SAXException
	        {
	        if(name.equals("DNA"))
	            {
	            this.baos=new ByteArrayOutputStream(this.reserve);
	            }
	        }
		
	    @Override
	    public void characters(char[] ch, int start, int length)
	            throws SAXException
	        {
	        if(this.baos==null) return;
	        for(int i=0;i< length;++i)
	            {
	            char c= Character.toUpperCase(ch[start+i]);
	            if(Character.isWhitespace(c)) continue;
	            this.baos.write((byte)c);
	            }
	        }
	
	
		public GenomicSequence getSequence(String chrom, int chromStart0, int chromEnd0)
				throws IOException
			{
			if(chromStart0 <0 || chromStart0 >=chromEnd0)
				{
				throw new IllegalArgumentException("Error in start/end");
				}
			this.reserve=(1+(chromEnd0-chromStart0));
			this.baos=null;
			try
				{
				String uri="http://genome.ucsc.edu/cgi-bin/das/"+
						BackLocate.this.genomeVersion+
						"/dna?segment="+
						URLEncoder.encode(chrom+":"+(chromStart0+1)+","+(chromEnd0+2), "UTF-8")
						;
				info(uri);
				InputStream in=new URL(uri).openStream();
				this.parser.parse(in, this);
				in.close();
				GenomicSequence g= new GenomicSequence(
					this.baos.toByteArray(),
					chrom,
					chromStart0
					);
				this.baos=null;
				return g;
				}
			catch (SAXException err)
				{
				throw new IOException(err);
				}
			
			}
		
		}

	
	

	

		

	private void backLocate(
		KnownGene gene,
		String geneName,
		char aa1,char aa2,
		int peptidePos1
		) throws IOException
		{
		final int extra=1000;
		
		GeneticCode geneticCode=getGeneticCodeByChromosome(gene.getChromosome());
		RNASequence wildRNA=null;
		ProteinCharSequence wildProt=null;
		
	        		
	        		
		if(genomicSeq==null ||
	               !gene.getChromosome().equals(genomicSeq.getChrom()) ||
	               !(genomicSeq.getChromStart()<=gene.getTxStart() && gene.getTxEnd()<= genomicSeq.getChromEnd())
	               )
        	{
        	this.info("fetch genome");
        	this.genomicSeq=this.dasServer.getSequence(
    	            		gene.getChromosome(),
    	            		Math.max(gene.getTxStart()-extra,0),
    	            		gene.getTxEnd()+extra
    	            		);
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
	
	private void run(LineIterator in) throws IOException,SQLException
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
			Set<String> kgIds=new HashSet<String>();
			PreparedStatement pstmt=connection.prepareStatement("select kgID from kgXref where geneSymbol=?");
			pstmt.setString(1, geneName);
			ResultSet row=pstmt.executeQuery();
			while(row.next())
				{
				kgIds.add(row.getString(1));
				}
			row.close();
			pstmt.close();
			if(kgIds.isEmpty())
				{
				warning("No kgXref found for "+geneName);
				continue;
				}
			StringBuilder sqlquery=new StringBuilder("select ");
			for(int i=0;i< KnownGene.SQL_COLUMNS.length;++i)
				{
				if(i>0) sqlquery.append(",");
				sqlquery.append(KnownGene.SQL_COLUMNS[i]);
				}
			sqlquery.append(" from knownGene where name=?");
			
			pstmt=connection.prepareStatement(sqlquery.toString());
			for(String kgId:kgIds)
				{
				pstmt.setString(1, kgId);
				row=pstmt.executeQuery();
				while(row.next())
					{
					KnownGene kg=new KnownGene(row);
					backLocate(kg, geneName, aa1, aa2, position1);
					}
				row.close();
				}
			pstmt.close();
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
	
	@Override
	public void printOptions(PrintStream out)
		{
		System.out.println(" -b ucsc.build default:"+ this.genomeVersion);
		System.out.println(" -p print mRNA & protein sequences");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		try {			
			com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
			int c;
			while((c=opt.getopt(args,getGetOptDefault()+ "b:p"))!=-1)
				{
				switch(c)
					{
					case 'b': this.genomeVersion=opt.getOptArg();break;
					case 'p': this.printSequences=true;break;
					default: 
						{
						switch(handleOtherOptions(c, opt, null))
							{
							case EXIT_FAILURE: return -1;
							case EXIT_SUCCESS: return 0;
							default: break;
							}
						}
					}
				}

			
			try
				{
				Class.forName("com.mysql.jdbc.Driver");
				}
			catch(Exception err)
				{
				error(err);
				error(getMessageBundle("cannot.load.mysql.driver"));
				return -1;
				}
			this.connection=DriverManager.getConnection(
				"jdbc:mysql://genome-mysql.cse.ucsc.edu/"+
				this.genomeVersion+
				"?user=genome&password="
				);
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
			CloserUtil.close(connection);
			}	
		}
	public static void main(String[] args)
		{
		new BackLocate().instanceMainWithExit(args);
		}
	}
