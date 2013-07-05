/**
 * Author:
 * 	Pierre Lindenbaum PhD
 * Date:
 * 	Dec-2010
 * Contact:
 * 	plindenbaum@yahoo.fr
 * Reference:
 *   http://plindenbaum.blogspot.com/2011/01/my-tool-to-annotate-vcf-files.html
 * Motivation:
 * 	Annotate a VCF file with the UCSC data. No SQL Driver required or a local database.
 * Compilation:
 */
package com.github.lindenb.jvarkit.tools.vcfannot;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.Closeable;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.StringReader;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLEncoder;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.util.Log;

import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.LineReader;
import org.broad.tribble.readers.TabixReader;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

import com.github.lindenb.jvarkit.util.AbstractVCFFilter;
import com.github.lindenb.jvarkit.util.GeneticCode;
import com.github.lindenb.jvarkit.util.picard.IOUtils;



/** a pair chromosome / position */
class ChromPosition
	implements Comparable<ChromPosition>
	{
	private String chromosome;
	private int position;
	public ChromPosition(String chromosome,int position)
		{
		this.chromosome=chromosome;
		this.position=position;
		}
	
	public String getChromosome()
		{
		return chromosome;
		}
	
	public int getPosition()
		{
		return position;
		}
	
	@Override
	public boolean equals(Object obj)
		{
		if (this == obj) { return true; }
		if (obj == null) { return false; }
		if (getClass() != obj.getClass()) { return false; }
		ChromPosition other = (ChromPosition) obj;
		if (position != other.position) { return false; }
		if (!chromosome.equalsIgnoreCase(other.chromosome)) { return false; }
		return true;
		}
	
	@Override
	public int hashCode()
		{
		final int prime = 31;
		int result = 1;
		result = prime * result + chromosome.hashCode();
		result = prime * result + position;
		return result;
		}
	
	@Override
	public int compareTo(ChromPosition o)
		{
		int i= getChromosome().compareToIgnoreCase(o.getChromosome());
		if(i!=0) return i;
		return getPosition()-o.getPosition();
		}
	@Override
	public String toString()
		{
		return getChromosome()+":"+getPosition();
		}
	}

/**
 * A record in a VCF file
 */
class VCFCall
	implements Comparable<VCFCall> //order by chromPosition
	{
	@SuppressWarnings("unused")
	private static Logger LOG=Logger.getLogger("vcf.annotator");
	/** the position */
	private ChromPosition chromPosition;
	/** columns in the VCF */
	private String columns[];
	
	/** cstor */
	VCFCall(String columns[])
		{
		this.columns=columns;
		if(!columns[0].toLowerCase().startsWith("chr"))
			{
			columns[0]="chr"+columns[0];
			}
		this.chromPosition=new ChromPosition(columns[0], Integer.parseInt(columns[1]));
		}
	
	/** get the columns from the VCF line */
	public String[] getColumns()
		{
		return columns;
		}
	
	/** get the position */
	public ChromPosition getChromPosition()
		{
		return this.chromPosition;
		}
	/** compare by position */
	@Override
	public int compareTo(VCFCall o)
		{
		return getChromPosition().compareTo(o.getChromPosition());
		}
	/** returns the VCF line */
	public String getLine()
		{
		StringBuilder line=new StringBuilder(100);
		for(int i=0;i< columns.length;++i)
			{
			if(i!=0) line.append("\t");
			line.append(columns[i]);
			}
		return line.toString();
		}
	
	
	@Override
	public String toString()
		{
		return getLine();
		}
	
	public void addProperties(String opcode,Map<String, String> map)
		{
		StringBuilder b=new StringBuilder();
		boolean first=true;
		for(String key:map.keySet())
			{
			if(!first) b.append("|");
			first=false;
			b.append(key);
			b.append(":");
			b.append(map.get(key));
			}
		addProperty(opcode, b.toString());
		}
	
	public void addProperty(String key,String value)
		{
		if(columns[7].equals(".")) columns[7]="";	
		if(!columns[7].isEmpty()) this.columns[7]+=";";
		columns[7]+=(key+"="+value);
		}
	
	public void addId(String newId)
		{
		String rsId= this.getColumns()[2];
		if(rsId.equals(".")) rsId="";
		Set<String> set=new HashSet<String>(Arrays.asList(rsId.split("[;]")));
		
		set.remove("");
		set.add(newId);
		rsId="";
		for(String s:set)
			{
			if(!rsId.isEmpty()) rsId+=";";
			rsId+=s;
			}
		this.getColumns()[2]=rsId;
		}
	}

/**
 * A VCF file
 * @author pierre
 *
 */
class VCFFile
	{
	private static Logger LOG=Logger.getLogger("vcf.annotator");
	//private static final String DEFAULT_HEADER="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample";
	private List<String> headers=new ArrayList<String>();
	private List<VCFCall> calls=new ArrayList<VCFCall>(10000);
	
	public VCFFile()
		{
		
		}
	
	public List<String> getHeaders()
		{
		return headers;
		}
	
	public List<VCFCall> getCalls()
		{
		return calls;
		}
	
	public  int lowerBound( ChromPosition position)
		{
		return lowerBound(0, getCalls().size(), position);
		}
    /** C+ lower_bound */
    public  int lowerBound(
                int first, int last,
                ChromPosition position
                )
        {
        int len = last - first;
        while (len > 0)
                {
                int half = len / 2;
                int middle = first + half;
                VCFCall call= getCalls().get(middle);
                if ( call.getChromPosition().compareTo(position) < 0  )
                        {
                        first = middle + 1;
                        len -= half + 1;
                        }
                else
                        {
                        len = half;
                        }
                }
   
        return first;
        }
    /** get the calls at given position */
	public List<VCFCall> get(ChromPosition pos)
    	{
    	int i=lowerBound(0,getCalls().size(),pos);
    	List<VCFCall> array=new ArrayList<VCFCall>(5);
    	while(i< getCalls().size())
    		{
    		VCFCall call=getCalls().get(i);
    		if(!call.getChromPosition().equals(pos)) break;
    		array.add(call);
    		++i;
    		}
    	return array;
    	}
	
	/** prints the VCF line */
	public void print(PrintWriter out)
		{
		for(String header:getHeaders())
			{
			out.println(header);
			}
		for(VCFCall c:getCalls())
			{
			out.println(c.getLine());
			}
		out.flush();
		}
	
	/** read VCF file */
	private void read(BufferedReader in)
	throws IOException
		{
		Pattern tab=Pattern.compile("[\t]");
		String line;
		while((line=in.readLine())!=null)
			{
			LOG.info(line);
			if(!line.startsWith("#")) break;
			this.headers.add(line);
			}
		
		if(!headers.isEmpty())
			{
			final String fmt="##fileformat=VCFv";
			String first= headers.get(0);
			
			if(first.startsWith("##format"))
				{
				first="##file"+first.substring(2);
				}
			
			if(!(first.startsWith(fmt)))
				{
				throw new IOException("firt line should starts with "+fmt);
				}
			String last=headers.get(headers.size()-1);
			if(!last.startsWith("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"))
				{
				throw new IOException("Error in header got "+line+" but expected "+
						"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
				}
			}
		else
			{
			this.headers.add("##fileformat=VCFv4.0");
			this.headers.add("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample");
			}
		
		while(line!=null)
			{
			//LOG.info(line);
			if(line.startsWith("#")) throw new IOException("line starting with # after header!"+line);
			String tokens[]=tab.split(line);
			if(tokens.length<8) throw new IOException("illegal number of columns in "+line);
			getCalls().add(new VCFCall(tokens));
			line=in.readLine();
			}
		
		Collections.sort(getCalls());
		LOG.info("vcf:"+getCalls().size()+" calls");
		}
	
	public Set<String> getChromosomes()
		{
		Set<String> set=new HashSet<String>();
		for(VCFCall c: getCalls())
			{
			set.add(c.getChromPosition().getChromosome());
			}
		LOG.info(set.toString());
		return set;
		}
	
	public void addHeader(String key,String value)
		{
		while(!key.startsWith("##")) key="#"+key;
		String line= key+"="+value;
		if(this.headers.contains(line)) return;
		this.headers.add(this.headers.size()-1,line);
		}
	
	public void addInfo(
		String id,Integer number,String type,String desc)
		{
		String line;
		if(this.headers.get(0).startsWith("##fileformat=VCFv4"))
			{
			line="<ID="+id+",Number="+(number==null?".":number.toString())+
				",Type="+type+",Description=\""+desc+"\">";
			}
		else if(this.headers.get(0).startsWith("##fileformat=VCFv3"))
			{
			line=id+","+(number==null?".":number.toString())+
			","+type+",\""+desc+"\"";
			}
		else
			{
			throw new IllegalArgumentException("VCF format not handled. "+this.headers.get(0));
			}
		addHeader("##INFO",line);
		}
	
	
	
	public static VCFFile parse(BufferedReader in)
	throws IOException
		{
		VCFFile vcf=new VCFFile();
		vcf.read(in);
		return vcf;
		}
	}


/** CharSeq a simple string impl */
interface CharSeq
	{
	public int length();
	public char charAt(int i);
	}

/**
 * Abstract implementation of CharSeq
 */
abstract class AbstractCharSeq implements CharSeq
	{
	AbstractCharSeq()
		{
		}
	
	@Override
	public int hashCode()
		{
		return getString().hashCode();
		}
	
	public String getString()
		{
		StringBuilder b=new StringBuilder(length());
		for(int i=0;i< length();++i) b.append(charAt(i));
		return b.toString();
		}
	
	@Override
	public String toString()
		{
		return getString();
		}
	}




/**
 * 
 * A GenomicSequence
 *
 */
class GenomicSequence
	extends AbstractCharSeq
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


class DefaultCharSeq extends AbstractCharSeq
	{
	private CharSequence seq;
	DefaultCharSeq(CharSequence seq)
		{
		this.seq=seq;
		}
	@Override
	public char charAt(int i)
		{
		return seq.charAt(i);
		}
	@Override
	public int length() {
		return seq.length();
		}
	}

class MutedSequence extends AbstractCharSeq
	{
	private CharSequence wild;
	private Map<Integer, Character> pos2char=new TreeMap<Integer, Character>();
	MutedSequence(CharSequence wild)
		{
		this.wild=wild;
		}
	
	void put(int pos,char c)
		{
		this.pos2char.put(pos, c);
		}
	
	@Override
	public char charAt(int i)
		{
		Character c= pos2char.get(i);
		return c==null?wild.charAt(i):c;
		}
	
	@Override
	public int length()
		{
		return this.wild.length();
		}
	}


class ProteinCharSequence extends AbstractCharSeq
	{
	private CharSeq cDNA;
	private GeneticCode geneticCode;
	ProteinCharSequence(GeneticCode geneticCode,CharSeq cDNA)
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
 * 
 * KnownGene
 *
 */
class KnownGene
	{
	private String name;
	private String chrom;
	private char strand;
	private int txStart;
	private int txEnd;
	private int cdsStart;
	private int cdsEnd;
	private int exonStarts[];
	private int exonEnds[];
	private String geneSymbol;
	
	abstract class Segment
		{
		private int index;
		protected Segment(int index)
			{
			this.index=index;
			}
		
		public int getIndex()
			{
			return index;
			}
		
		public KnownGene getGene()
			{
			return KnownGene.this;
			}
		
		public boolean contains(int position)
			{
			return getStart()<=position && position< getEnd();
			}
		public abstract boolean isSplicingAcceptor(int position);
		public abstract boolean isSplicingDonor(int position);
		public boolean isSplicing(int position)
			{
			return isSplicingAcceptor(position) || isSplicingDonor(position);
			}
		
		public abstract String getName();
		public abstract int getStart();
		public abstract int getEnd();
		}
	
	class Exon extends Segment
		{
		private Exon(int index)
			{
			super(index);
			}
		
		@Override
		public String getName()
			{
			if(getGene().getStrand()=='+')
				{
				return "Exon "+(getIndex()+1);
				}
			else
				{
				return "Exon "+(getGene().getExonCount()-getIndex());
				}
			}
		
		@Override
		public int getStart()
			{
			return getGene().getExonStart(getIndex());
			}
		
		@Override
		public int getEnd()
			{
			return getGene().getExonEnd(getIndex());
			}
		
		@Override
		public String toString()
			{
			return getName();
			}
		
		
		public Intron getNextIntron()
			{
			if(getIndex()+1>=getGene().getExonCount()) return null;
			return getGene().getIntron(getIndex());
			}
		public Intron getPrevIntron()
			{
			if(getIndex()<=0) return null;
			return getGene().getIntron(getIndex()-1);
			}
		
		@Override
		public boolean isSplicingAcceptor(int position)
			{
			if(!contains(position)) return false;
			if(isForward())
				{
				if(getIndex()== 0) return false;
				return position==getStart();
				}
			else
				{
				if(getIndex()+1== getGene().getExonCount()) return false;
				return position==getEnd()-1;
				}
			}
		
		@Override
		public boolean isSplicingDonor(int position)
			{
			if(!contains(position)) return false;
			if(isForward())
				{
				if(getIndex()+1== getGene().getExonCount()) return false;
				return  (position==getEnd()-1) ||
						(position==getEnd()-2) ||
						(position==getEnd()-3) ;
				}
			else
				{
				if(getIndex()== 0) return false;
				return  (position==getStart()+0) ||
						(position==getStart()+1) ||
						(position==getStart()+2) ;
				}
			}
		
		}
		
	class Intron extends Segment
			{
			Intron(int index)
				{
				super(index);
				}
			
			@Override
			public int getStart()
				{
				return getGene().getExonEnd(getIndex());
				}
			
			@Override
			public int getEnd()
				{
				return getGene().getExonStart(getIndex()+1);
				}
			
			@Override
			public String getName() {
				if(getGene().isForward())
					{
					return "Intron "+(getIndex()+1);
					}
				else
					{
					return "Intron "+(getGene().getExonCount()-getIndex());
					}
				}

			public boolean isSplicingAcceptor(int position)
				{
				if(!contains(position)) return false;
				if(isForward())
					{
					return  (position==getEnd()-1) ||
							(position==getEnd()-2);
					}
				else
					{
					return	position==getStart() ||
							position==getStart()+1;
					}
				}
			

			public boolean isSplicingDonor(int position)
				{
				if(!contains(position)) return false;
				if(isForward())
					{
					return	position==getStart() ||
							position==getStart()+1;
							
					}
				else
					{
					return  (position==getEnd()-1) ||
							(position==getEnd()-2);
					}
				}
			
			}
	
		/**
		 * 
		 * KnownGene 
		 * 
		 */
		public KnownGene(String tokens[])
			throws IOException
			{
			this.name = tokens[0];
			this.geneSymbol=tokens[0];
			this.chrom= tokens[1];
	        this.strand = tokens[2].charAt(0);
	        this.txStart = Integer.parseInt(tokens[3]);
	        this.txEnd = Integer.parseInt(tokens[4]);
	        this.cdsStart= Integer.parseInt(tokens[5]);
	        this.cdsEnd= Integer.parseInt(tokens[6]);
	        int exonCount=Integer.parseInt(tokens[7]);
	        this.exonStarts = new int[exonCount];
	        this.exonEnds = new int[exonCount];
	            
            
            int index=0;
            for(String s: tokens[8].split("[,]"))
            	{
            	this.exonStarts[index++]=Integer.parseInt(s);
            	}
            index=0;
            for(String s: tokens[9].split("[,]"))
            	{
            	this.exonEnds[index++]=Integer.parseInt(s);
            	}
			}
		
		/** returns knownGene ID */
		public String getName()
			{
			return this.name;
			}
		
		/** returns chromosome name */
		public String getChromosome()
			{
			return this.chrom;
			}
		
		/** returns the strand */
		public char getStrand()
			{
			return strand;
			}
		boolean isForward()
        	{
        	return getStrand()=='+';
        	}

		public int getTxStart()
			{
			return this.txStart;
			}

		public int getTxEnd()
			{
			return this.txEnd;
			}
		

		public int getCdsStart()
			{
			return this.cdsStart;
			}
		

		public int getCdsEnd()
			{
			return this.cdsEnd;
			}
		

		public int getExonStart(int index)
			{
			return this.exonStarts[index];
			}
		

		public int getExonEnd(int index)
			{
			return this.exonEnds[index];
			}
		

		public Exon getExon(int index)
			{
			return new Exon(index);
			}
		public Intron getIntron(int i)
			{
			return new Intron(i);
			}
		public int getExonCount()
			{
			return this.exonStarts.length;
			}
		public String getGeneSymbol()
			{
			return geneSymbol;
			}
		
		public void setGeneSymbol(String geneSymbol)
			{
			this.geneSymbol = geneSymbol;
			}
		
		}


/*************************************************************************************/

/**
 * PredictionAnnotator
 *
 */
class PredictionAnnotator implements Closeable
	{
	static final String KEY_TYPE="type";
	static final String KEY_SPLICING="splicing";
	static Logger LOG=Logger.getLogger("vcf.annotator");
	private Map<String, List<KnownGene>> chrom2genes=new HashMap<String, List<KnownGene>>();
	private IndexedFastaSequenceFile indexedFastaSequenceFile;
	
	
	private static char complement(char c)
		{
		switch(c)
			{
			case 'A': return 'T';
			case 'T': return 'A';
			case 'G': return 'C';
			case 'C': return 'G';
			default:throw new IllegalArgumentException(""+c);
			}
		}
	
	

	TabixReader tabixReader=null;
	
	private List<KnownGene> getGenes(final VariantContext ctx)
	throws IOException	
		{
		Pattern tab=Pattern.compile("[\t]");
		List<KnownGene> L=new ArrayList<KnownGene>();
		String line;
		TabixReader.Iterator iter=this.tabixReader.query("");//TODO
		while((iter!=null && (line=iter.next())!=null))
			{
			L.add(new KnownGene(tab.split(line)));
			}
		return L;
		}
	
	@Override
	public void close()
		{
		if(this.tabixReader!=null) this.tabixReader.close();
		this.tabixReader=null;
		}
	
	
	public VariantContext run(VariantContext ctx) throws IOException
			{
			VariantContextBuilder b=new VariantContextBuilder(ctx);
			
            List<KnownGene> genes= getGenes(ctx);
            for(Allele a: ctx.getAlleles())
            	{
            	
            	}
            return ctx;
			}
	
	private VariantContext run(
			VariantContext ctx,Allele allele,
			final List<KnownGene> genes)
			throws IOException
			{
			String genomicSeq=null;
           int position= ctx.getStart()-1;

            String ref=ctx.getReference().getDisplayString();

			String alt=allele.getDisplayString();
        	
            if(ref.equals("A"))
    			{
    				 if(alt.equals("W")) { alt="T"; }
    			else if(alt.equals("M")) { alt="C"; }
    			else if(alt.equals("R")) { alt="G"; }
    			}
    		else if(ref.equals("C"))
    			{
    				 if(alt.equals("S")) { alt="G"; }
    			else if(alt.equals("M")) { alt="A"; }
    			else if(alt.equals("Y")) { alt="T"; }
    			}
    		else if(ref.equals("G"))
    			{
    				 if(alt.equals("S")) { alt="C"; }
    			else if(alt.equals("K")) { alt="T"; }
    			else if(alt.equals("R")) { alt="A"; }
    			}
    		else if(ref.equals("T"))
    			{
    				 if(alt.equals("W")) { alt="A"; }
    			else if(alt.equals("K")) { alt="G"; }
    			else if(alt.equals("Y")) { alt="C"; }
    			}
    		
    		LOG.info(ref+" "+alt);
            
            
            if(genes.isEmpty())
            	{
            	LOG.info("GENOMIC");
            	continue;
            	}
            
            for(KnownGene gene:genes)
            	{
            	LOG.info(gene.getName());
            	
            	//switch to 0 based coordinate
        		
            	Map<String, String> annotations=new HashMap<String, String>();
        		
        		
            	if( (ref.equals("A") || ref.equals("T") || ref.equals("G") || ref.equals("C")) &&
            		(alt.equals("A") || alt.equals("T") || alt.equals("G") || alt.equals("C"))
            		)
	        		{
	        		LOG.info("fetch genome");
            		GeneticCode geneticCode=GeneticCode.getByChromosome(gene.getChromosome());
	        		StringBuilder wildRNA=null;
	        		ProteinCharSequence wildProt=null;
	        		ProteinCharSequence mutProt=null;
	        		MutedSequence mutRNA=null;
	        		int position_in_cdna=-1;
	        		
	        		
	        		if(genomicSeq==null ||
	        	               !gene.getChromosome().equals(genomicSeq.getChrom()) ||
	        	               !(genomicSeq.getChromStart()<=gene.getTxStart() && gene.getTxEnd()<= genomicSeq.getChromEnd())
	        	               )
    	            	{
    	            	final int maxTry=20;
    	            	for(int tryCount=1;tryCount<=maxTry;++tryCount)
    	            		{
    	            		genomicSeq=null;
    	            		try
    	            			{
    	            			genomicSeq=this.indexedFastaSequenceFile.getSequence(
			    	            		gene.getChromosome(),
			    	            		Math.max(gene.getTxStart()+1,0),
			    	            		gene.getTxEnd()
			    	            		);
	    	          
    	            			}
    	            		catch (Exception e)
    	            			{
								LOG.info("Cannot get DAS-DNA sequence "+e.getMessage());
								genomicSeq=null;
								}
    	            		if(genomicSeq!=null)
    	            			{
    	            			break;
    	            			}
    	            		LOG.info("try to get DAS-DNA "+(tryCount)+"/"+maxTry);
    	            		}
    	            	if(genomicSeq==null)
    	            		{
    	            		throw new IOException("Cannot get DAS-DNA");
    	            		}
    	            	}
	        		
	        		if(!String.valueOf(genomicSeq.charAt(position)).equalsIgnoreCase(ref))
	        			{
	        			System.err.println("Warning REF!=GENOMIC SEQ!!! at "+genomicSeq.charAt(position)+"/"+ref);
	        			return;
	        			}
	        		
	        		if(gene.isForward())
	            		{
	            		if(position < gene.getCdsStart())
	            			{
	            			annotations.put("type", "UTR5");
	            			}
	            		else if( gene.getCdsEnd()<= position )
	            			{
	            			annotations.put("type", "UTR3");
	            			}
	            		
	            		int exon_index=0;
	            		while(exon_index< gene.getExonCount())
	            			{
	            			KnownGene.Exon exon= gene.getExon(exon_index);
	            			for(int i= exon.getStart();
	            					i< exon.getEnd();
	            					++i)
	            				{
	            				if(i==position)
	        						{
	        						annotations.put("exon.name", exon.getName());
	        						}
	            				if(i< gene.getCdsStart()) continue;
	            				if(i>=gene.getCdsEnd()) break;
	        					
	        					if(wildRNA==null)
	        						{
	        						wildRNA=new StringBuilder();
	        						mutRNA=new MutedSequence(wildRNA);
	        						}
	        					
	        					if(i==position)
	        						{
	        						annotations.put("type", "EXON");
	        						annotations.put("exon.name",exon.getName());
	        						position_in_cdna=wildRNA.length();
	        						annotations.put("pos.cdna", String.valueOf(position_in_cdna));
	        						//in splicing ?
	        						if(exon.isSplicing(position))
	        							{
	        							annotations.put(KEY_SPLICING, "SPLICING");
	        							
	        							if(exon.isSplicingAcceptor(position))
	        								{
	        								annotations.put(KEY_SPLICING, "SPLICING_ACCEPTOR");
	        								}
	        							else  if(exon.isSplicingDonor(position))
	        								{
	        								annotations.put(KEY_SPLICING, "SPLICING_DONOR");
	        								}
	        							}
	        						}
	        					
	            				wildRNA.append(genomicSeq.charAt(i));
	            				
	            				if(i==position)
	            					{
	            					mutRNA.put(
	            							position_in_cdna,
	            							alt.charAt(0)
	            							);
	            					}
	            				
	            				if(wildRNA.length()%3==0 && wildRNA.length()>0 && wildProt==null)
		            				{
		            				wildProt=new ProteinCharSequence(geneticCode,new DefaultCharSeq(wildRNA));
		            				mutProt=new ProteinCharSequence(geneticCode,mutRNA);
		            				}
	            				}
	            			KnownGene.Intron intron= exon.getNextIntron();
	            			if(intron!=null && intron.contains(position))
	            				{
	            				annotations.put("intron.name",intron.getName());
	            				annotations.put("type", "INTRON");
	            				
	            				if(intron.isSplicing(position))
	        						{
	            					annotations.put(KEY_SPLICING, "INTRON_SPLICING");
	        						if(intron.isSplicingAcceptor(position))
	        							{
	        							annotations.put(KEY_SPLICING, "INTRON_SPLICING_ACCEPTOR");
	        							}
	        						else if(intron.isSplicingDonor(position))
	        							{
	        							annotations.put(KEY_SPLICING, "INTRON_SPLICING_DONOR");
	        							}
	        						}
	            				}
	            			++exon_index;
	            			}
	            		
	            		
	            		
	            		}
	            	else // reverse orientation
	            		{
	            	
	            		if(position < gene.getCdsStart())
	            			{
	            			annotations.put(KEY_TYPE, "UTR3");
	            			}
	            		else if( gene.getCdsEnd()<=position )
	            			{
	            			annotations.put(KEY_TYPE, "UTR5");
	            			}
	            	
	            		
	            		int exon_index = gene.getExonCount()-1;
	            		while(exon_index >=0)
	            			{
	            			KnownGene.Exon exon= gene.getExon(exon_index);
	            			for(int i= exon.getEnd()-1;
	            				    i>= exon.getStart();
	            				--i)
	            				{
	            				if(i==position)
	        						{
	            					annotations.put("exon.name", exon.getName());
	        						}
	            				if(i>= gene.getCdsEnd()) continue;
	            				if(i<  gene.getCdsStart()) break;
	            				
	            				if(wildRNA==null)
	        						{
	        						wildRNA=new StringBuilder();
	        						mutRNA=new MutedSequence(wildRNA);
	        						}
	            				
	            				if(i==position)
	        						{
	            					annotations.put(KEY_TYPE, "EXON");
	            					position_in_cdna=wildRNA.length();
	        						annotations.put("pos.cdna",String.valueOf(position_in_cdna));
	        						//in splicing ?
	        						if(exon.isSplicing(position))
	        							{
	        							annotations.put(KEY_SPLICING, "INTRON_SPLICING");
	        							
	        							if(exon.isSplicingAcceptor(position))
	        								{
	        								annotations.put(KEY_SPLICING, "INTRON_SPLICING_ACCEPTOR");
	        								}
	        							else  if(exon.isSplicingDonor(position))
	        								{
	        								annotations.put(KEY_SPLICING, "INTRON_SPLICING_DONOR");
	        								}
	        							}
	        						
	        						
	        						mutRNA.put(
	        								position_in_cdna,
	        								complement(alt.charAt(0))
	        								);
	        						}
	            				
	            				wildRNA.append(complement(genomicSeq.charAt(i)));
	            				if( wildRNA.length()%3==0 &&
	            					wildRNA.length()>0 &&
	            					wildProt==null)
		            				{
		            				wildProt=new ProteinCharSequence(geneticCode,new DefaultCharSeq(wildRNA));
		            				mutProt=new ProteinCharSequence(geneticCode,mutRNA);
		            				}
	            				
	            				}
	            			
	            			KnownGene.Intron intron= exon.getPrevIntron();
	            			if(intron!=null &&
	            				intron.contains(position))
	            				{
	            				annotations.put("intron.name",intron.getName());
	            				annotations.put(KEY_TYPE, "INTRON");
	            				
	            				if(intron.isSplicing(position))
	        						{
	            					annotations.put(KEY_SPLICING, "INTRON_SPLICING");
	        						if(intron.isSplicingAcceptor(position))
	        							{
	        							annotations.put(KEY_SPLICING, "INTRON_SPLICING_ACCEPTOR");
	        							}
	        						else if(intron.isSplicingDonor(position))
	        							{
	        							annotations.put(KEY_SPLICING, "INTRON_SPLICING_DONOR");
	        							}
	        						}
	            				}
	            			--exon_index;
	            			}

	            		}//end of if reverse
	        		if( wildProt!=null &&
	        			mutProt!=null && 
	        			position_in_cdna>=0)
		    			{
	            		int pos_aa=position_in_cdna/3;
	            		int mod= position_in_cdna%3;
	            		annotations.put("wild.codon",""+
	            			wildRNA.charAt(position_in_cdna-mod+0)+
	            			wildRNA.charAt(position_in_cdna-mod+1)+
	            			wildRNA.charAt(position_in_cdna-mod+2)
	            			);
	            		annotations.put("mut.codon",""+
	            			mutRNA.charAt(position_in_cdna-mod+0)+
	            			mutRNA.charAt(position_in_cdna-mod+1)+
	            			mutRNA.charAt(position_in_cdna-mod+2)
	            			);
	            		annotations.put("position.protein",String.valueOf(pos_aa+1));
	            		annotations.put("wild.aa",String.valueOf(wildProt.charAt(pos_aa)));
	            		annotations.put("mut.aa",String.valueOf(mutProt.charAt(pos_aa)));
		    			if(isStop(wildProt.charAt(pos_aa)) &&
		    			   !isStop(mutProt.charAt(pos_aa)))
		    				{
		    				annotations.put("type", "EXON_STOP_LOST");
		    				}
		    			else if( !isStop(wildProt.charAt(pos_aa)) &&
		    				 isStop(mutProt.charAt(pos_aa)))
		    				{
		    				annotations.put("type", "EXON_STOP_GAINED");
		    				}
		    			else if(wildProt.charAt(pos_aa)==mutProt.charAt(pos_aa))
		    				{
		    				annotations.put("type", "EXON_CODING_SYNONYMOUS");
		    				}
		    			else
		    				{
		    				annotations.put("type", "EXON_CODING_NON_SYNONYMOUS");
		    				}
		    			}
	        		}//end of simpe ATCG
            	else
            		{
	        		Integer wildrna=null;
	        		int position_in_cdna=-1;
	        		
	        		
	        		
	        		if(gene.isForward())
	            		{
	            		if(position < gene.getCdsStart())
	            			{
	            			annotations.put("type", "UTR5");
	            			}
	            		else if( gene.getCdsEnd()<= position )
	            			{
	            			annotations.put("type", "UTR3");
	            			}
	            		
	            		int exon_index=0;
	            		while(exon_index< gene.getExonCount())
	            			{
	            			KnownGene.Exon exon= gene.getExon(exon_index);
	            			for(int i= exon.getStart();
	            					i< exon.getEnd();
	            					++i)
	            				{
	            				if(i==position)
	        						{
	        						annotations.put("exon.name", exon.getName());
	        						}
	            				if(i< gene.getCdsStart()) continue;
	            				if(i>=gene.getCdsEnd()) break;
	        					
	        					if(wildrna==null)
	        						{
	        						wildrna=0;
	        						}
	        					
	        					if(i==position)
	        						{
	        						annotations.put(KEY_TYPE, "EXON");
	        						annotations.put("exon.name",exon.getName());
	        						position_in_cdna=wildrna;
	        						annotations.put("pos.cdna", String.valueOf(wildrna));
	        						//in splicing ?
	        						if(exon.isSplicing(position))
	        							{
	        							annotations.put(KEY_SPLICING, "SPLICING");
	        							
	        							if(exon.isSplicingAcceptor(position))
	        								{
	        								annotations.put(KEY_SPLICING, "SPLICING_ACCEPTOR");
	        								}
	        							else  if(exon.isSplicingDonor(position))
	        								{
	        								annotations.put(KEY_SPLICING, "SPLICING_DONOR");
	        								}
	        							}
	        						}
	        					
	        					wildrna++;
	            				}
	            			KnownGene.Intron intron= exon.getNextIntron();
	            			if(intron!=null && intron.contains(position))
	            				{
	            				annotations.put("intron.name",intron.getName());
	            				annotations.put("type", "INTRON");
	            				
	            				if(intron.isSplicing(position))
	        						{
	            					annotations.put(KEY_SPLICING, "INTRON_SPLICING");
	        						if(intron.isSplicingAcceptor(position))
	        							{
	        							annotations.put(KEY_SPLICING, "INTRON_SPLICING_ACCEPTOR");
	        							}
	        						else if(intron.isSplicingDonor(position))
	        							{
	        							annotations.put(KEY_SPLICING, "INTRON_SPLICING_DONOR");
	        							}
	        						}
	            				}
	            			++exon_index;
	            			}
	            		}
	            	else // reverse orientation
	            		{
	            	
	            		if(position < gene.getCdsStart())
	            			{
	            			annotations.put(KEY_TYPE, "UTR3");
	            			}
	            		else if( gene.getCdsEnd()<=position )
	            			{
	            			annotations.put(KEY_TYPE, "UTR5");
	            			}
	            	
	            		
	            		int exon_index = gene.getExonCount()-1;
	            		while(exon_index >=0)
	            			{
	            			KnownGene.Exon exon= gene.getExon(exon_index);
	            			for(int i= exon.getEnd()-1;
	            				    i>= exon.getStart();
	            				--i)
	            				{
	            				if(i==position)
	        						{
	            					annotations.put("exon.name", exon.getName());
	        						}
	            				if(i>= gene.getCdsEnd()) continue;
	            				if(i<  gene.getCdsStart()) break;
	            				
	            				if(wildrna==null)
	        						{
	        						wildrna=0;
	        						}
	            				
	            				if(i==position)
	        						{
	            					annotations.put(KEY_TYPE, "EXON");
	            					position_in_cdna=wildrna;
	        						annotations.put("pos.cdna",String.valueOf(position_in_cdna));
	        						//in splicing ?
	        						if(exon.isSplicing(position))
	        							{
	        							annotations.put(KEY_SPLICING, "INTRON_SPLICING");
	        							
	        							if(exon.isSplicingAcceptor(position))
	        								{
	        								annotations.put(KEY_SPLICING, "INTRON_SPLICING_ACCEPTOR");
	        								}
	        							else  if(exon.isSplicingDonor(position))
	        								{
	        								annotations.put(KEY_SPLICING, "INTRON_SPLICING_DONOR");
	        								}
	        							}
	        						}
	            				
	            				wildrna++;
	            				}
	            			
	            			KnownGene.Intron intron= exon.getPrevIntron();
	            			if(intron!=null &&
	            				intron.contains(position))
	            				{
	            				annotations.put("intron.name",intron.getName());
	            				annotations.put(KEY_TYPE, "INTRON");
	            				
	            				if(intron.isSplicing(position))
	        						{
	            					annotations.put(KEY_SPLICING, "INTRON_SPLICING");
	        						if(intron.isSplicingAcceptor(position))
	        							{
	        							annotations.put(KEY_SPLICING, "INTRON_SPLICING_ACCEPTOR");
	        							}
	        						else if(intron.isSplicingDonor(position))
	        							{
	        							annotations.put(KEY_SPLICING, "INTRON_SPLICING_DONOR");
	        							}
	        						}
	            				}
	            			--exon_index;
	            			}

	            		}//end of if reverse
	        		if( wildrna!=null &&
	        			position_in_cdna>=0)
		    			{
	            		int pos_aa=position_in_cdna/3;
	            		annotations.put("position.protein",String.valueOf(pos_aa+1));
		    			}
            		}
            	annotations.put("strand", ""+gene.getStrand());
            	annotations.put("kgId", gene.getName());
            	annotations.put("geneSymbol", gene.getGeneSymbol());
            	call.addProperties("PREDICTION", annotations);
            	}
            
	    return b.make();
		}
	
	
	private boolean isStop(char aa)
		{
		return aa=='*';
		}
	
	
	}

/**
 * VCFAnnotator
 * Annotator for VCF
 *
 */
public class VCFAnnotator extends AbstractVCFFilter
	{
	static Log LOG=Log.getInstance(VCFAnnotator.class);
	public File IN;
	
	public VCFAnnotator()
		{
		}
	
	@Override
	protected void doWork(LineReader in, VariantContextWriter out)
		throws IOException
		{
		// TODO Auto-generated method stub
		
		VCFCodec vcfCodec=new VCFCodec();
		vcfCodec.readHeader(in);
		String line;
		while((line=in.readLine())!=null)
			{
			VariantContext ctx=vcfCodec.decode(line);
			
			}
			
		}
	
	public static void main(String[] args)
		{
		new VCFAnnotator().instanceMainWithExit(args);
		}
	}