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


History:
* 2015 creation copied from http://plindenbaum.blogspot.fr/2010/11/blastxmlannotations.html

*/
package com.github.lindenb.jvarkit.tools.biostar;

/*
BEGIN_DOC

## Example

```
$ cat ~/jeter.blastn.xml 
<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
(...)
<Hit>
  <Hit_num>1</Hit_num>
  <Hit_id>gi|14971104|gb|AF338247.1|</Hit_id>
  <Hit_def>Human rotavirus A strain M clone M1 NSP3 genes, complete cds</Hit_def>
  <Hit_accession>AF338247</Hit_accession>
  <Hit_len>2032</Hit_len>
  <Hit_hsps>
    <Hsp>
      <Hsp_num>1</Hsp_num>
```

```
$ java -jar dist/biostar3654.jar ~/jeter.blastn.xml 2> /dev/null  | cut -c-${COLUMNS} 

QUERY: No definition line
       ID:Query_186611 Len:980
>Human rotavirus A strain M clone M1 NSP3 genes, complete cds
 AF338247
 id:gi|14971104|gb|AF338247.1| len:2032

   e-value:0 gap:0 bitScore:1764.98

QUERY 000000001 GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGC 000000050
                ||||||||||||||||||||||||||||||||||||||||||||||||||
HIT   000000001 GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGC 000000050
                ################################################## source organi
                ##################################                 5'UTR
                                                  ################ CDS codon_sta

QUERY 000000051 AGATGGTAAGCTCTATTATTAATACTTCTTTTGAAGCTGCAGTCGTTGCT 000000100
                ||||||||||||||||||||||||||||||||||||||||||||||||||
HIT   000000051 AGATGGTAAGCTCTATTATTAATACTTCTTTTGAAGCTGCAGTCGTTGCT 000000100
                ################################################## source organi
                ################################################## CDS codon_sta
(...)

```

END_DOC
*/
import gov.nih.nlm.ncbi.blast.Hit;
import gov.nih.nlm.ncbi.blast.Hsp;
import gov.nih.nlm.ncbi.blast.Iteration;
import gov.nih.nlm.ncbi.insdseq.INSDFeature;
import gov.nih.nlm.ncbi.insdseq.INSDFeatureIntervals;
import gov.nih.nlm.ncbi.insdseq.INSDInterval;
import gov.nih.nlm.ncbi.insdseq.INSDQualifier;
import htsjdk.samtools.util.CloserUtil;

import java.util.List;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.net.URL;
import java.util.ArrayList;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Unmarshaller;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLResolver;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.XMLEvent;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.ncbi.NcbiApiKey;
import com.github.lindenb.jvarkit.util.ncbi.NcbiConstants;
/**
BEGIN_DOC

## Example

```
$ cat ~/jeter.blastn.xml 
<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
(...)
<Hit>
  <Hit_num>1</Hit_num>
  <Hit_id>gi|14971104|gb|AF338247.1|</Hit_id>
  <Hit_def>Human rotavirus A strain M clone M1 NSP3 genes, complete cds</Hit_def>
  <Hit_accession>AF338247</Hit_accession>
  <Hit_len>2032</Hit_len>
  <Hit_hsps>
    <Hsp>
      <Hsp_num>1</Hsp_num>
```

```
$ java -jar dist/biostar3654.jar ~/jeter.blastn.xml 2> /dev/null  | cut -c-${COLUMNS} 

QUERY: No definition line
       ID:Query_186611 Len:980
>Human rotavirus A strain M clone M1 NSP3 genes, complete cds
 AF338247
 id:gi|14971104|gb|AF338247.1| len:2032

   e-value:0 gap:0 bitScore:1764.98

QUERY 000000001 GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGC 000000050
                ||||||||||||||||||||||||||||||||||||||||||||||||||
HIT   000000001 GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGC 000000050
                ################################################## source organi
                ##################################                 5'UTR
                                                  ################ CDS codon_sta

QUERY 000000051 AGATGGTAAGCTCTATTATTAATACTTCTTTTGAAGCTGCAGTCGTTGCT 000000100
                ||||||||||||||||||||||||||||||||||||||||||||||||||
HIT   000000051 AGATGGTAAGCTCTATTATTAATACTTCTTTTGAAGCTGCAGTCGTTGCT 000000100
                ################################################## source organi
                ################################################## CDS codon_sta
(...)
```

END_DOC

 */

@Program(name="biostar3654",
description="show blast alignment with annotations",
	biostars=3654,
	keywords={"blast","xml","annotation"})
public class Biostar3654 extends Launcher
	{
	private static final Logger LOG=Logger.build(Biostar3654.class).make();
	
	@SuppressWarnings("unused")
	private static final gov.nih.nlm.ncbi.blast.ObjectFactory _fool_javac1=null;
	@SuppressWarnings("unused")
	private static final gov.nih.nlm.ncbi.insdseq.ObjectFactory _fool_javac2=null;
	
	
	private PrintWriter pw=null;
	/** left margin */
	private int margin=9;
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile=null;
	/** length of a fasta line */
	@Parameter(names={"-L","--length"},description="Fasta Line kength")
	private int fastaLineLength=50;
	@ParametersDelegate
	private NcbiApiKey ncbiApiKey = new NcbiApiKey();

	
	/** XML input factory */
	private XMLInputFactory xif;
	/** transforms XML/DOM to GBC entry */
	private Unmarshaller unmarshaller;
	
	
	
	
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
			for(final INSDFeature feature:this.features)
				{
				if(feature.getINSDFeatureIntervals()==null) continue;
				if(feature.getINSDFeatureIntervals().getINSDInterval()==null) continue;
				//loop over the coordinates
				for(final INSDInterval interval:feature.getINSDFeatureIntervals().getINSDInterval())
					{
					int intervalFrom=0;
					int intervalTo=0;
					//is it an interval ?
					if( interval.getINSDIntervalFrom()!=null &&
						interval.getINSDIntervalTo()!=null )
						{
						intervalFrom = Integer.parseInt(interval.getINSDIntervalFrom());
						intervalTo   = Integer.parseInt(interval.getINSDIntervalTo());
						}
					//is it a single point ?
					else if(interval.getINSDIntervalPoint()!=null &&
						(intervalFrom=Integer.parseInt(interval.getINSDIntervalPoint()))>=this.seqStart &&
						intervalFrom< this.seqEnd
						)
						{
						intervalFrom = Integer.parseInt(interval.getINSDIntervalPoint());
						intervalTo   = intervalFrom;
						}
					else
						{
						continue;
						}
					if(intervalFrom> intervalTo)
						{
						int tmp=intervalFrom;
						intervalFrom=intervalTo;
						intervalTo=tmp;
						}
					intervalTo++;
					
					
					if(intervalFrom>this.seqEnd) continue;
					if(intervalTo<this.seqStart) continue;
					
					
					//margin left
					Biostar3654.this.pw.printf("      %"+margin+"s ","");
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
							Biostar3654.this.pw.print(isSeq?(isGap?":":"#"):"-");
							}
						else //not in the feature
							{
							Biostar3654.this.pw.print(" ");
							}
						//extends the current position if current char is a base/aminoacid
						if(Character.isLetter(getSequence().charAt(this.stringStart+i)))
							{
							genome+=this.sign;
							}
						}
					Biostar3654.this.pw.print(" ");
					Biostar3654.this.pw.print(feature.getINSDFeatureKey());
					//Biostar3654.this.pw.print(" ");
					//Biostar3654.this.pw.print(feature.getINSDFeatureLocation());//no because using seq_start & seq_stop with efetch change this
					//print the infos
					if( feature.getINSDFeatureQuals()!=null && feature.getINSDFeatureQuals().getINSDQualifier()!=null ) {
						for(final INSDQualifier qual:feature.getINSDFeatureQuals().getINSDQualifier())
							{
							Biostar3654.this.pw.print(" ");
							Biostar3654.this.pw.print(qual.getINSDQualifierName());
							Biostar3654.this.pw.print(":");
							Biostar3654.this.pw.print(qual.getINSDQualifierValue());
							}
						}
					
					Biostar3654.this.pw.println();
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
	private List<INSDFeature> fetchAnnotations(
		final String database,
		final String acn,int start,int end)
		throws Exception
		{
		InputStream in=null;
		XMLEventReader r=null;
		final List<INSDFeature> L=new ArrayList<INSDFeature>();
		if(start>end) return fetchAnnotations(database,acn,end,start);
		try 
			{
			
			if(acn!=null && !acn.isEmpty() && !acn.startsWith("Query"))
				{
				String uri=
						NcbiConstants.efetch()+"?db="+database+
						"&id="+StringUtils.escapeHttp(acn)+
						"&rettype=gbc&retmode=xml&seq_start="+start+"&seq_stop="+end+
						this.ncbiApiKey.getAmpParamValue()
						;
				LOG.info(uri);
				in = new URL(uri).openStream();
				r=this.xif.createXMLEventReader(in);
				while(r.hasNext())
					{
					XMLEvent  evt=r.peek();
					if(evt.isStartElement() &&
							evt.asStartElement().getName().getLocalPart().equals("INSDFeature"))
							{
							INSDFeature feature = this.unmarshaller.unmarshal(r,
									INSDFeature.class).getValue();
							INSDFeatureIntervals its = feature.getINSDFeatureIntervals();
							if(its==null || its.getINSDInterval().isEmpty()) continue;
							for(INSDInterval interval:its.getINSDInterval())
								{
								//when using seq_start and seq_stop , the NCBI shifts the data...
								if( interval.getINSDIntervalFrom()!=null &&
									interval.getINSDIntervalTo()!=null)
									{
									interval.setINSDIntervalFrom(String.valueOf(Integer.parseInt(interval.getINSDIntervalFrom())+start-1));
									interval.setINSDIntervalTo(String.valueOf(Integer.parseInt(interval.getINSDIntervalTo())+start-1));
									}
								else if( interval.getINSDIntervalPoint()!=null)
									{
									interval.setINSDIntervalPoint(String.valueOf(Integer.parseInt(interval.getINSDIntervalPoint())+start-1));
									}
								}
							L.add(feature);
							}
						else
							{
							r.next();//consumme
							}			
					}
				
				}
			}
		catch(Exception err) {
			LOG.error(err);
			}
		finally
			{
			CloserUtil.close(r);
			CloserUtil.close(in);
			}
			LOG.info("N(INSDFeature)="+L.size());
			//not found, return empty table
			return L;
			}
	
	/** parses BLAST output */
	private void parseBlast(XMLEventReader r)  throws Exception
		{
		String database="nucleotide";
		while(r.hasNext())
			{
			XMLEvent evt=r.peek();
			if(evt.isStartElement() &&
					evt.asStartElement().getName().getLocalPart().equals("BlastOutput_program"))
				{
				r.next();
				String BlastOutput_program = r.getElementText();
				if("blastn".equals(BlastOutput_program))
					{
					database="nucleotide";
					}
				else if("blastp".equals(BlastOutput_program) )
					{
					database="protein";
					}
				else
					{
					throw new IOException("only blastn && blastn are supported: "+database);
					}
				}
			else if(evt.isStartElement() &&
				evt.asStartElement().getName().getLocalPart().equals("Iteration"))
				{
				Iteration iteration= this.unmarshaller.unmarshal(r, Iteration.class).getValue();
				parseIteration(database,iteration);
				}
			else
				{
				r.next();//consumme
				}
			}
		}
	
	private void parseIteration(String database,Iteration iteration) throws Exception
		{
		this.pw.println("QUERY: "+iteration.getIterationQueryDef());
		this.pw.println("       ID:"+iteration.getIterationQueryID()+" Len:"+iteration.getIterationQueryLen());
		
		for(Hit hit:iteration.getIterationHits().getHit())
			{
			
			this.pw.println(">"+hit.getHitDef());
			this.pw.println(" "+hit.getHitAccession());
			this.pw.println(" id:"+hit.getHitId()+" len:"+hit.getHitLen());
			for(Hsp hsp :hit.getHitHsps().getHsp())
				{
				List<INSDFeature> qFeatures= fetchAnnotations(
						database,
						iteration.getIterationQueryID(),
						Integer.parseInt(hsp.getHspQueryFrom()),
						Integer.parseInt(hsp.getHspQueryTo())
						);
				List<INSDFeature> hFeatures= fetchAnnotations(
						database,
						hit.getHitAccession(),
						Integer.parseInt(hsp.getHspHitFrom()),
						Integer.parseInt(hsp.getHspHitTo())
						);

				this.pw.println();
				this.pw.println("   e-value:"+hsp.getHspEvalue()+" gap:"+hsp.getHspGaps()+" bitScore:"+hsp.getHspBitScore());
				this.pw.println();
				//create the Printer for the Query and the Hit
				QPrinter qPrinter=new QPrinter(hsp,qFeatures);
				HPrinter hPrinter=new HPrinter(hsp,hFeatures);
				
				//loop over the lines
				while(qPrinter.next() && hPrinter.next())
					{
					qPrinter.print();
					this.pw.printf("QUERY %0"+margin+"d ",qPrinter.seqStart);
					this.pw.print(hsp.getHspQseq().substring(qPrinter.stringStart,qPrinter.stringEnd));
					this.pw.printf(" %0"+margin+"d",qPrinter.seqEnd-(qPrinter.sign));
					this.pw.println();
					this.pw.printf("      %"+margin+"s ","");
					this.pw.print(hsp.getHspMidline().substring(qPrinter.stringStart,qPrinter.stringEnd));
					this.pw.println();
					this.pw.printf("HIT   %0"+margin+"d ",hPrinter.seqStart);
					this.pw.print(hsp.getHspHseq().substring(hPrinter.stringStart,hPrinter.stringEnd));
					this.pw.printf(" %0"+margin+"d",hPrinter.seqEnd-(hPrinter.sign));
					this.pw.println();
					hPrinter.print();
					this.pw.println();
					}
				this.pw.flush();
				}
			
			
		}
		
		
		//System.err.println("OK");
		}

	@Override
	public int doWork(final List<String> args) {
		
		if(!this.ncbiApiKey.isApiKeyDefined()) {
			LOG.error("NCBI API key is not defined");
			return -1;
			}
		
		try
			{
			
			//create a Unmarshaller for genbank
			JAXBContext jc = JAXBContext.newInstance(
					"gov.nih.nlm.ncbi.insdseq:gov.nih.nlm.ncbi.blast");
			this.unmarshaller=jc.createUnmarshaller();
	
			this.xif=XMLInputFactory.newFactory();
			this.xif.setXMLResolver(new XMLResolver()
				{
				@Override
				public Object resolveEntity(String publicID,
						String systemID, String baseURI, String namespace)
						throws XMLStreamException {
							return new ByteArrayInputStream(new byte[0]);
						}
				});
			this.pw = super.openFileOrStdoutAsPrintWriter(this.outputFile);
			
			//read from stdin
			if(args.isEmpty())
				{
				XMLEventReader r=this.xif.createXMLEventReader(stdin(), "UTF-8");
				this.parseBlast(r);
				r.close();
				}
			else
				{
				//loop over the files
				for(String inputName:args)
					{
					LOG.info("Reading "+inputName);
					FileReader fr=new java.io.FileReader(inputName);
					XMLEventReader r=this.xif.createXMLEventReader(fr);
					this.parseBlast(r);
					r.close();
					fr.close();
					}
				}
			pw.flush();
			return 0;
			}
		catch(Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.pw);
			}
		}
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new Biostar3654().instanceMainWithExit(args);
		}

	}
