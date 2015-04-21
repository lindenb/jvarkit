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
* 2015 creation copied from http://plindenbaum.blogspot.fr/2010/11/blastxmlannotations.html

*/
package com.github.lindenb.jvarkit.tools.biostar;

import gov.nih.nlm.ncbi.blast.BlastOutput;
import gov.nih.nlm.ncbi.blast.BlastOutputIterations;
import gov.nih.nlm.ncbi.blast.Hit;
import gov.nih.nlm.ncbi.blast.Hsp;
import gov.nih.nlm.ncbi.blast.Iteration;
import gov.nih.nlm.ncbi.insdseq.INSDFeature;
import gov.nih.nlm.ncbi.insdseq.INSDInterval;
import gov.nih.nlm.ncbi.insdseq.INSDQualifier;
import gov.nih.nlm.ncbi.insdseq.INSDSeq;
import gov.nih.nlm.ncbi.insdseq.INSDSeqFeatureTable;
import gov.nih.nlm.ncbi.insdseq.INSDSet;

import java.io.PrintStream;
import java.util.List;
import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Unmarshaller;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.xml.sax.EntityResolver;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

import com.github.lindenb.jvarkit.knime.AbstractKnimeApplication;

public class Biostar3654 extends AbstractKnimeApplication
	{
	@SuppressWarnings("unused")
	private static final gov.nih.nlm.ncbi.blast.ObjectFactory _fool_javac1=null;
	@SuppressWarnings("unused")
	private static final gov.nih.nlm.ncbi.insdseq.ObjectFactory _fool_javac2=null;
	
	/** left margin */
	private int margin=9;
	/** length of a fasta line */
	private int fastaLineLength=50;
	/** xml parser */
	private DocumentBuilder docBuilder;
	/** transforms XML/DOM to GBC entry */
	private Unmarshaller gbcUnmarshaller;
	
	/** abstract class writing the line of an alignment */
	private abstract class AbstractHspPrinter
		{
		/** index in sequence for this line*/
		int seqStart;
		/** end of sequence for this line */
		int seqEnd;
		/** forward/reverse */
		int sign;
		/** index in alignment string */
		int stringStart;
		/** end index in alignment string */
		int stringEnd;
		/** current HSP */
		protected Hsp hsp;
		/** features to be printed */
		private List<INSDFeature> features;
		
		/** get sequence to be displayed */
		public abstract String getSequence();
		/** get starting index of the sequence */
		public abstract int getSeqFrom();
		/** get end index of the sequence */
		public abstract int getSeqTo();
		
		protected AbstractHspPrinter( Hsp hsp, List<INSDFeature> features)
			{
			this.hsp=hsp;
			this.features=features;
			
			this.sign= (getSeqFrom()< getSeqTo()?1:-1);
			this.seqStart=getSeqFrom();
			this.seqEnd=this.seqStart;
			this.stringStart=0;
			this.stringEnd=0;
			}
		/** can we print another line ? yes ? init the data */
		boolean next()
			{
			if(this.stringEnd>=getSequence().length()) return false;
			this.seqStart=this.seqEnd;
			this.stringStart=this.stringEnd;
			for(int i=0;i< fastaLineLength &&
					this.stringStart+i< getSequence().length();
					++i)
				{
				if(Character.isLetter(getSequence().charAt(this.stringStart+i)))
					{
					this.seqEnd+=this.sign;
					}
				this.stringEnd++;
				}
			return true;
			}
		
		/** print the  line */
		void print()
			{
			/* loop over the feature */
			for(INSDFeature feature:this.features)
				{
				if(feature.getINSDFeatureIntervals()==null) continue;
				//loop over the coordinates
				for(INSDInterval interval:feature.getINSDFeatureIntervals().getINSDInterval())
					{
					int intervalFrom=0;
					int intervalTo=0;
					//is it an interval ?
					if( interval.getINSDIntervalFrom()!=null &&
						interval.getINSDIntervalTo()!=null &&
						!(
						  (intervalFrom=Integer.parseInt(interval.getINSDIntervalFrom()))>=this.seqEnd ||
						  (intervalTo=Integer.parseInt(interval.getINSDIntervalTo()))<this.seqStart
						))
						{
						intervalTo++;
						}
					//is it a single point ?
					else if(interval.getINSDIntervalPoint()!=null &&
						(intervalFrom=Integer.parseInt(interval.getINSDIntervalPoint()))>=this.seqStart &&
						intervalFrom< this.seqEnd
						)
						{
						intervalTo=intervalFrom+1;
						}
					else
						{
						continue;
						}
					if(intervalFrom> intervalTo)
						{
						//uhh ???
						continue;
						}
					//margin left
					System.out.printf("      %"+margin+"s ","");
					//trim the positions
					intervalFrom=Math.max(this.seqStart,intervalFrom);
					intervalTo=Math.min(this.seqEnd,intervalTo);
					int genome=this.seqStart;
					
					//loop over the line
					for(	int i=0;i< fastaLineLength &&
							this.stringStart+i< this.stringEnd;
							++i)
						{
						boolean isSeq=Character.isLetter(getSequence().charAt(this.stringStart+i));
						boolean isGap=hsp.getHspMidline().charAt(this.stringStart+i)==' ';
						//in the feature
						if(intervalFrom<=genome && genome< intervalTo)
							{
							System.out.print(isSeq?(isGap?":":"#"):"-");
							}
						else //not in the feature
							{
							System.out.print(" ");
							}
						//extends the current position if current char is a base/aminoacid
						if(Character.isLetter(getSequence().charAt(this.stringStart+i)))
							{
							genome+=this.sign;
							}
						}
					System.out.print(" ");
					System.out.print(feature.getINSDFeatureKey());
					System.out.print(" ");
					System.out.print(feature.getINSDFeatureLocation());
					//print the infos
					for(INSDQualifier qual:feature.getINSDFeatureQuals().getINSDQualifier())
						{
						System.out.print(" ");
						System.out.print(qual.getINSDQualifierName());
						System.out.print(":");
						System.out.print(qual.getINSDQualifierValue());
						}
					
					System.out.println();
					}
				}
			}
		}
	
	/** specialized AbstractHspPrinter for the QUERY */
	private	class QPrinter
		extends AbstractHspPrinter
		{
		QPrinter( Hsp hsp, List<INSDFeature> features)
			{
			super(hsp,features);
			}
		@Override
		public int getSeqFrom()
			{
			return Integer.parseInt(this.hsp.getHspQueryFrom());
			}
		@Override
		public int getSeqTo()
			{
			return Integer.parseInt(this.hsp.getHspQueryTo());
			}
		public String getSequence()
			{
			return this.hsp.getHspQseq();
			}
		}
	
	/** specialized AbstractHspPrinter for the HIT */
	private	class HPrinter
	extends AbstractHspPrinter
		{
		HPrinter( Hsp hsp, List<INSDFeature> features)
			{
			super(hsp,features);
			}
		
		@Override
		public int getSeqFrom()
			{
			return Integer.parseInt(this.hsp.getHspHitFrom());
			}
		
		@Override
		public int getSeqTo()
			{
			return Integer.parseInt(this.hsp.getHspHitTo());
			}
		
		public String getSequence()
			{
			return this.hsp.getHspHseq();
			}
		}
	
	
	
	/** fetches the annotation for a given entry if the name starts with gi|.... */
	private List<INSDFeature> fetchAnnotations(String name)
		throws Exception
		{
		int pipe;
		if(name!=null &&
		name.startsWith("gi|") &&
		(pipe=name.indexOf('|',3))!=-1)
			{
			String gi=name.substring(3,pipe);
			String uri="http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id="+gi+"&rettype=gbc&retmode=xml";
			INSDSet set=INSDSet.class.cast(this.gbcUnmarshaller.unmarshal(this.docBuilder.parse(uri)));
			if(!set.getINSDSeq().isEmpty())
				{
				INSDSeq seq= set.getINSDSeq().get(0);
				INSDSeqFeatureTable table=seq.getINSDSeqFeatureTable();
				return table.getINSDFeature();
				}
			}
		//not found, return empty table
		return new ArrayList<INSDFeature>();
		}
	
	/** parses BLAST output */
	private void parseBlast(BlastOutput blast)  throws Exception
		{
		System.out.println("QUERY: "+blast.getBlastOutputQueryDef());
		System.out.println("       ID:"+blast.getBlastOutputQueryID()+" Len:"+blast.getBlastOutputQueryLen());
		List<INSDFeature> qFeatures= fetchAnnotations(blast.getBlastOutputQueryID());
		BlastOutputIterations iterations=blast.getBlastOutputIterations();
		for(Iteration iteration:iterations.getIteration())
			{
			for(Hit hit:iteration.getIterationHits().getHit())
				{
				System.out.println(">"+hit.getHitDef());
				System.out.println(" "+hit.getHitAccession());
				System.out.println(" id:"+hit.getHitId()+" len:"+hit.getHitLen());
				List<INSDFeature> hFeatures= fetchAnnotations(hit.getHitId());
				for(Hsp hsp :hit.getHitHsps().getHsp())
					{
					System.out.println();
					System.out.println("   e-value:"+hsp.getHspEvalue()+" gap:"+hsp.getHspGaps()+" bitScore:"+hsp.getHspBitScore());
					System.out.println();
					//create the Printer for the Query and the Hit
					QPrinter qPrinter=new QPrinter(hsp,qFeatures);
					HPrinter hPrinter=new HPrinter(hsp,hFeatures);
					
					//loop over the lines
					while(qPrinter.next() && hPrinter.next())
						{
						qPrinter.print();
						System.out.printf("QUERY %0"+margin+"d ",qPrinter.seqStart);
						System.out.print(hsp.getHspQseq().substring(qPrinter.stringStart,qPrinter.stringEnd));
						System.out.printf(" %0"+margin+"d",qPrinter.seqEnd-(qPrinter.sign));
						System.out.println();
						System.out.printf("      %"+margin+"s ","");
						System.out.print(hsp.getHspMidline().substring(qPrinter.stringStart,qPrinter.stringEnd));
						System.out.println();
						System.out.printf("HIT   %0"+margin+"d ",hPrinter.seqStart);
						System.out.print(hsp.getHspHseq().substring(hPrinter.stringStart,hPrinter.stringEnd));
						System.out.printf(" %0"+margin+"d",hPrinter.seqEnd-(hPrinter.sign));
						System.out.println();
						hPrinter.print();
						System.out.println();
						}
					}
				
				}
			}
		
		
		//System.err.println("OK");
		}

	@Override
	public int executeKnime(List<String> args)
		{
		try
			{
			//create a DOM parser
			DocumentBuilderFactory f=DocumentBuilderFactory.newInstance();
			f.setCoalescing(true);
			f.setNamespaceAware(true);
			f.setValidating(false);
			f.setExpandEntityReferences(true);
			f.setIgnoringComments(false);
			f.setIgnoringElementContentWhitespace(true);
			this.docBuilder= f.newDocumentBuilder();
			this.docBuilder.setEntityResolver(new EntityResolver()
				{
				@Override
				public InputSource resolveEntity(String publicId, String systemId)
						throws SAXException, IOException
					{
					return new InputSource(new StringReader(""));
					}
				});
			//create a Unmarshaller for genbank
			JAXBContext jc = JAXBContext.newInstance("gov.nih.nlm.ncbi.insdseq");
			this.gbcUnmarshaller=jc.createUnmarshaller();
	
			
			jc = JAXBContext.newInstance("gov.nih.nlm.ncbi.blast");
			Unmarshaller unmarshaller=jc.createUnmarshaller();
			//read from stdin
			if(args.isEmpty())
				{
				this.parseBlast(BlastOutput.class.cast(unmarshaller.unmarshal(this.docBuilder.parse(System.in))));
				}
			else
				{
				//loop over the files
				for(String inputName:args)
					{
					this.parseBlast(BlastOutput.class.cast(unmarshaller.unmarshal(this.docBuilder.parse(new File(inputName)))));
					}
				}
			return 0;
			}
		catch(Throwable err)
			{
			error(err);
			return -1;
			}
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println("-o (fileout). Default: stdout");
		out.println("-L (int) line length: default: "+this.fastaLineLength);
		super.printOptions(out);
		}

	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"o:L:"))!=-1)
			{
			switch(c)
				{
				case 'L': this.fastaLineLength=Math.max(1,Integer.parseInt(opt.getOptArg()));break;
				case 'o': setOutputFile(opt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}

		return mainWork(opt.getOptInd(), args);
		}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new Biostar3654().instanceMainWithExit(args);
		}

	}
