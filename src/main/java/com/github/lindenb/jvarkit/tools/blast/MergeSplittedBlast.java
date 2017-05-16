package com.github.lindenb.jvarkit.tools.blast;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import gov.nih.nlm.ncbi.blast.Hit;
import gov.nih.nlm.ncbi.blast.HitHsps;
import gov.nih.nlm.ncbi.blast.Hsp;
import gov.nih.nlm.ncbi.blast.Iteration;
import gov.nih.nlm.ncbi.blast.IterationHits;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLEventWriter;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLResolver;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.XMLEvent;

import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.blast.BlastHspAlignment;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
/**
 * 
BEGIN_DOC
 

## Example

### Creating the blast database
the long sequence must be splitted using a sliding+overlapping window.
The fasta  header must be formatted as:
```
(chromName):(1-based-start-position)-(1-based-end-position):(chromLength)
```
### Makefile for test:
```make
BLASTBIN=/commun/data/packages/ncbi/ncbi-blast-2.2.28+/bin

all:blast.xml

blast.xml:query.fa database.fa
        ${BLASTBIN}/blastn -db database.fa -query query.fa -outfmt 5 -out $@ -dust no

query.fa:
        echo -e ">q1\nCCTACGTTCAATATTACAGGCGAACATACCTACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCGCTTTCCACACAGACATCATAACAAA" > $@

database.fa:
        curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chrM.fa.gz" |\
        gunzip -c | grep -v '^>' | tr -d "\n" |\
        awk '{L=length($$0);for(i=1;i+77<=L;i+=65) printf(">chrM:%d-%d:%d\n%s\n",i,i+77,L,substr($$0,i,77));}' |\
        fold -w 50 > $@ && \
        ${BLASTBIN}/makeblastdb -in $@ -dbtype nucl
```
The blast.xml produced is:
```xml
<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_version>BLASTN 2.2.28+</BlastOutput_version>
  <BlastOutput_reference>Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb Miller (2000), "A greedy algorithm for aligning DNA sequences", J Comput Biol 2000; 7(1-2):203-14.</BlastOutput_reference>
  <BlastOutput_db>database.fa</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>q1</BlastOutput_query-def>
  <BlastOutput_query-len>123</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_expect>10</Parameters_expect>
      <Parameters_sc-match>1</Parameters_sc-match>
      <Parameters_sc-mismatch>-2</Parameters_sc-mismatch>
      <Parameters_gap-open>0</Parameters_gap-open>
      <Parameters_gap-extend>0</Parameters_gap-extend>
      <Parameters_filter>m;</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>q1</Iteration_query-def>
      <Iteration_query-len>123</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|3</Hit_id>
          <Hit_def>chrM:196-273:16571</Hit_def>
          <Hit_accession>3</Hit_accession>
          <Hit_len>77</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>143.312</Hsp_bit-score>
              <Hsp_score>77</Hsp_score>
              <Hsp_evalue>1.29152e-37</Hsp_evalue>
              <Hsp_query-from>31</Hsp_query-from>
              <Hsp_query-to>107</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>77</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>77</Hsp_identity>
              <Hsp_positive>77</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>77</Hsp_align-len>
              <Hsp_qseq>TACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCGCTTTCCACA</Hsp_qseq>
              <Hsp_hseq>TACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCGCTTTCCACA</Hsp_hseq>
              <Hsp_midline>|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>2</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|2</Hit_id>
          <Hit_def>chrM:131-208:16571</Hit_def>
          <Hit_accession>2</Hit_accession>
          <Hit_len>77</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>78.6796</Hsp_bit-score>
              <Hsp_score>42</Hsp_score>
              <Hsp_evalue>3.69397e-18</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>42</Hsp_query-to>
              <Hsp_hit-from>36</Hsp_hit-from>
              <Hsp_hit-to>77</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>42</Hsp_identity>
              <Hsp_positive>42</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>42</Hsp_align-len>
              <Hsp_qseq>CCTACGTTCAATATTACAGGCGAACATACCTACTAAAGTGTG</Hsp_qseq>
              <Hsp_hseq>CCTACGTTCAATATTACAGGCGAACATACCTACTAAAGTGTG</Hsp_hseq>
              <Hsp_midline>||||||||||||||||||||||||||||||||||||||||||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>3</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|4</Hit_id>
          <Hit_def>chrM:261-338:16571</Hit_def>
          <Hit_accession>4</Hit_accession>
          <Hit_len>77</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>52.8265</Hsp_bit-score>
              <Hsp_score>28</Hsp_score>
              <Hsp_evalue>2.23898e-10</Hsp_evalue>
              <Hsp_query-from>96</Hsp_query-from>
              <Hsp_query-to>123</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>28</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>28</Hsp_identity>
              <Hsp_positive>28</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>28</Hsp_align-len>
              <Hsp_qseq>CCGCTTTCCACACAGACATCATAACAAA</Hsp_qseq>
              <Hsp_hseq>CCGCTTTCCACACAGACATCATAACAAA</Hsp_hseq>
              <Hsp_midline>||||||||||||||||||||||||||||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
      <Iteration_stat>
        <Statistics>
          <Statistics_db-num>254</Statistics_db-num>
          <Statistics_db-len>19558</Statistics_db-len>
          <Statistics_hsp-len>13</Statistics_hsp-len>
          <Statistics_eff-space>1788160</Statistics_eff-space>
          <Statistics_kappa>0.46</Statistics_kappa>
          <Statistics_lambda>1.28</Statistics_lambda>
          <Statistics_entropy>0.85</Statistics_entropy>
        </Statistics>
      </Iteration_stat>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
```
### Process the XML with MergeSplittedBlast:
```bash
$ java -jar dist/mergesplittedblast.jar blast.xml |\
  xmllint --format - 
```
### Output:
```xml
<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_version>BLASTN 2.2.28+</BlastOutput_version>
  <BlastOutput_reference>Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb Miller (2000), "A greedy algorithm for aligning DNA sequences", J Comput Biol 2000; 7(1-2):203-14.</BlastOutput_reference>
  <BlastOutput_db>database.fa</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>q1</BlastOutput_query-def>
  <BlastOutput_query-len>123</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_expect>10</Parameters_expect>
      <Parameters_sc-match>1</Parameters_sc-match>
      <Parameters_sc-mismatch>-2</Parameters_sc-mismatch>
      <Parameters_gap-open>0</Parameters_gap-open>
      <Parameters_gap-extend>0</Parameters_gap-extend>
      <Parameters_filter>m;</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>q1</Iteration_query-def>
      <Iteration_query-len>123</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|3</Hit_id>
          <Hit_def>chrM</Hit_def>
          <Hit_accession>3</Hit_accession>
          <Hit_len>16571</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>78.6796</Hsp_bit-score>
              <Hsp_score>123</Hsp_score>
              <Hsp_evalue>3.69397e-18</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>123</Hsp_query-to>
              <Hsp_hit-from>166</Hsp_hit-from>
              <Hsp_hit-to>288</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>42</Hsp_identity>
              <Hsp_positive>42</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>123</Hsp_align-len>
              <Hsp_qseq>CCTACGTTCAATATTACAGGCGAACATACCTACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCGCTTTCCACACAGACATCATAACAAA</Hsp_qseq>
              <Hsp_hseq>CCTACGTTCAATATTACAGGCGAACATACCTACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATAACAATTGAATGTCTGCACAGCCGCTTTCCACACAGACATCATAACAAA</Hsp_hseq>
              <Hsp_midline>|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
```
END_DOC
 */
@Program(name="mergesplittedblast",
	description="merge blast Hits from splitted BLAST database",
		biostars=90186,
		keywords={"blast"})
public class MergeSplittedBlast extends Launcher {
	private static final Logger LOG=Logger.build(MergeSplittedBlast.class).make();
	
	private Unmarshaller unmarshaller;
	private Marshaller marshaller;
	
	/* force javac to compile */
	@SuppressWarnings("unused")
	private gov.nih.nlm.ncbi.blast.ObjectFactory _ignore_for_javac=null;
	
	private static class Split
		{
		String chrom;
		int start1;
		int end1;
		int len;
		}
	
	
	
	private MergeSplittedBlast()
		{
		
		}
	
	private Split parse(Hit hit)
		{
		String s=hit.getHitDef();
		if(s==null ) return null;
		int colon=s.indexOf(':');
		if(colon==-1) return null;
		int dash=s.indexOf('-', colon+1);
		if(dash==-1) return null;
		int colon2=s.indexOf(':', dash+1);
		if(colon2==-1) return null;
		Split split=new Split();
		try	{
			split.chrom=s.substring(0, colon).trim();
			if(split.chrom.isEmpty()) return null;
			
			split.start1=Integer.parseInt(s.substring(colon+1,dash));
			if(split.start1<1) return null;
			split.end1=Integer.parseInt(s.substring(dash+1,colon2));
			if(split.end1<=split.start1) return null;
			split.len=Integer.parseInt(s.substring(colon2+1));
			}
		catch(NumberFormatException err)
			{
			LOG.error("expect (chrom:start-end:length) Cannot parse Hit.def: "+s ,err);
			return null;
			}
		return split;
		}
	
	
	private boolean overlap(int a_start1,int a_end1,int b_start1,int b_end1)
		{
		//debug(" overlapp ["+a_start1+"-"+a_end1+"] vs ["+b_start1+"-"+b_end1+"]");
		if(b_end1+1 < a_start1) return false;
		if(a_end1+1 < b_start1) return false;
		//debug(" ok overlap");
		return true;
		}
	private Hsp merge(Hsp hh0,Hsp hh1)
		{
		BlastHspAlignment aln0=new BlastHspAlignment(hh0);
		BlastHspAlignment aln1=new BlastHspAlignment(hh1);
		
		//order hsp on query pos
		if(aln0.getQueryFrom1() > aln1.getQueryFrom1())
			{
			Hsp th=hh0;
			hh0=hh1;
			hh1=th;
			aln0=new BlastHspAlignment(hh0);
			aln1=new BlastHspAlignment(hh1);
			}
		
		/* not the same strand */
		if(aln0.isPlusPlus()!=aln1.isPlusPlus())
			{
			//debug("strand");
			return null;
			}
		
		/* hits overlap */	
		if(!overlap(aln0.getHitFrom1(),aln0.getHitTo1(),aln1.getHitFrom1(),aln1.getHitTo1()))
			{
			//debug("no overlap hit");
			return null;
			}
		/* query overlap */
		if(!overlap(aln0.getQueryFrom1(),aln0.getQueryTo1(),aln1.getQueryFrom1(),aln1.getQueryTo1()))
			{
			//debug("no overlap query");
			return null;
			}
		//hit1 is contained in hit0
		if(aln0.getQueryFrom1() <= aln1.getQueryFrom1() && aln1.getQueryTo1()<=aln0.getQueryTo1())
			{
			//debug("contained");
			return hh0;
			}
		
		StringBuilder qsb=new StringBuilder();
		StringBuilder msb=new StringBuilder();
		StringBuilder hsb=new StringBuilder();
		int expect_hit=-1;
		int found_hit=-1;
		for(BlastHspAlignment.Align a:aln0)
			{
			if(a.getQueryIndex1()>=aln1.getQueryFrom1())
				{
				//debug("###BREAK###############");
				expect_hit=a.getHitIndex1();
				break;
				}
			qsb.append(a.getQueryChar());
			msb.append(a.getMidChar());
			hsb.append(a.getHitChar());
			}
		if(expect_hit==-1)
			{
			//debug("HU?");
			return null;
			}
		for(BlastHspAlignment.Align a:aln1)
			{
			if(a.getQueryIndex1()<aln1.getQueryFrom1()) continue;
			if(found_hit==-1)
				{
				found_hit=a.getHitIndex1();
				if(expect_hit!=found_hit)
					{
					LOG.info("Not the expected hit position "+expect_hit+"/"+found_hit);
					return null;
					}
				}
			qsb.append(a.getQueryChar());
			msb.append(a.getMidChar());
			hsb.append(a.getHitChar());
			}
		//info("\n"+qsb+"\n"+msb+"\n"+hsb);
		
		Hsp newHsp=BlastHspAlignment.cloneHsp(hh0);
		newHsp.setHspAlignLen(String.valueOf(msb.length()));
		newHsp.setHspMidline(msb.toString());
		newHsp.setHspQseq(qsb.toString());
		newHsp.setHspHseq(hsb.toString());
		newHsp.setHspQueryFrom(String.valueOf(Math.min(aln0.getQueryFrom1(),aln1.getQueryFrom1())));
		newHsp.setHspQueryTo(String.valueOf(Math.max(aln0.getQueryTo1(),aln1.getQueryTo1())));
		if(aln0.isPlusPlus())
			{
			newHsp.setHspHitFrom(String.valueOf(Math.min(aln0.getHitFrom1(),aln1.getHitFrom1())));
			newHsp.setHspHitTo(String.valueOf(Math.max(aln0.getHitTo1(),aln1.getHitTo1())));
			}
		else
			{
			newHsp.setHspHitFrom(String.valueOf(Math.max(aln0.getHitFrom1(),aln1.getHitFrom1())));
			newHsp.setHspHitTo(String.valueOf(Math.min(aln0.getHitTo1(),aln1.getHitTo1())));
			}
		newHsp.setHspGaps(String.valueOf(newHsp.getHspMidline().replaceAll("[^ ]", "").length()));
		newHsp.setHspScore(String.valueOf(newHsp.getHspMidline().replaceAll("[ ]", "").length()));
		//debug("success");
		return newHsp;
		}
	
	
	
	private Hit merge(List<Hit> hits)
		{
		Hit first=hits.get(0);
		Split firstSplit=parse(first);
		Hit newHit=new Hit();
		newHit.setHitAccession(first.getHitAccession());
		newHit.setHitDef(firstSplit.chrom);
		newHit.setHitId(first.getHitId());
		newHit.setHitNum(first.getHitNum());
		newHit.setHitLen(String.valueOf(firstSplit.len));
		newHit.setHitHsps(new HitHsps());
		
		List<Hsp> hsps=new ArrayList<Hsp>();
		
		
		for(Hit h0:hits)
			{
			/* fix hit_from/hit_to */
			Split split=parse(h0);
			List<Hsp> L=h0.getHitHsps().getHsp();
			for(int i=0;i< L.size();++i)
				{
				Hsp newHsp=BlastHspAlignment.cloneHsp(L.get(i));
				BlastHspAlignment aln=new BlastHspAlignment(newHsp);
				int h_from=aln.getHitFrom1();
				int h_to=aln.getHitTo1();
				
				h_from+=(split.start1-1);
				h_to+=(split.start1-1);
					
				
				newHsp.setHspHitFrom(String.valueOf(h_from));
				newHsp.setHspHitTo(String.valueOf(h_to));
				
				hsps.add(newHsp);
				}
			}
		
		
		boolean done=false;
		while(!done)
			{
			done=true;
			for(int i=0;i+1< hsps.size() ;++i)
				{
				Hsp hsp0=hsps.get(i);
				for(int j=i+1;j< hsps.size();++j)
					{
					//debug("comparing hsp "+i+" vs "+j+" N="+hsps.size());
					Hsp hsp1=hsps.get(j);
					Hsp newHitHsp=merge(hsp0,hsp1);
					if(newHitHsp!=null)
						{
						hsps.set(i,newHitHsp);
						hsps.remove(j);
						done=false;
						break;
						}
					}
				if(!done) break;
				}
			}
		for(int i=0;i< hsps.size();++i) hsps.get(i).setHspNum(String.valueOf(i+1));
		newHit.getHitHsps().getHsp().addAll(hsps);
		return newHit;
		}
	
	private Iteration merge(Iteration iteration)
		{
		if(iteration.getIterationHits().getHit().size()<=1) return iteration;
		Map<String,List<Hit>> chrom2hits=new LinkedHashMap<String,List<Hit>>();
		for(Hit hit:iteration.getIterationHits().getHit())
			{
			Split s=parse(hit);
			if(s==null) return iteration;
			List<Hit> L= chrom2hits.get(s.chrom);
			if(L==null)
				{
				L=new ArrayList<Hit>();
				 chrom2hits.put(s.chrom,L);
				}
			L.add(hit);
			}
		Iteration newiteration=new Iteration();
		List<Hit> newHits=new ArrayList<Hit>();
		for(String chrom:chrom2hits.keySet())
			{
			List<Hit> L= chrom2hits.get(chrom);
			Hit newHit=merge(L);
			newHits.add(newHit);
			}

		newiteration.setIterationIterNum(iteration.getIterationIterNum());
		newiteration.setIterationQueryID(iteration.getIterationQueryID());
		newiteration.setIterationQueryLen(iteration.getIterationQueryLen());
		newiteration.setIterationQueryDef(iteration.getIterationQueryDef());
		newiteration.setIterationMessage(iteration.getIterationMessage());
		newiteration.setIterationHits(new IterationHits());
		newiteration.getIterationHits().getHit().addAll(newHits);
		return newiteration;
		}
	
	private void run(
			XMLEventReader r,
			XMLEventWriter w
			)
			throws XMLStreamException,JAXBException
		{
		while(r.hasNext())
			{
			XMLEvent evt=r.peek();
			if(!(evt.isStartElement()))
				{
				w.add(r.nextEvent());
				continue;
				}
			String name=evt.asStartElement().getName().getLocalPart();
			if(name.equals("BlastOutput_program"))
				{
				w.add(r.nextEvent());
				evt=r.nextEvent();//TEXT
				if(!evt.isCharacters()) throw new XMLStreamException("Expected a string after "+name);
				String prog=evt.asCharacters().getData();
				if(!"blastn".equals(prog))
					{
					throw new XMLStreamException("Not a blastn input:"+prog);
					}
				w.add(evt);
				continue;
				}
			
			if(!name.equals("Iteration"))
				{
				w.add(r.nextEvent());
				continue;
				}
			
			Iteration iter= this.unmarshaller.unmarshal(r, Iteration.class).getValue();
			iter=merge(iter);
			this.marshaller.marshal(iter, w);
			}
		}
	
	@Override
	public int doWork(List<String> args) {
		XMLEventReader rx=null;
		XMLEventWriter wx=null;
		try
			{
			
			JAXBContext jc = JAXBContext.newInstance("gov.nih.nlm.ncbi.blast");
			this.unmarshaller=jc.createUnmarshaller();
			this.marshaller=jc.createMarshaller();
			this.marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT,true);
			this.marshaller.setProperty(Marshaller.JAXB_FRAGMENT,true);
			XMLInputFactory xmlInputFactory=XMLInputFactory.newFactory();
			xmlInputFactory.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, Boolean.FALSE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_COALESCING, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_REPLACING_ENTITY_REFERENCES, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_SUPPORTING_EXTERNAL_ENTITIES, Boolean.FALSE);
			xmlInputFactory.setXMLResolver(new XMLResolver()
				{
				@Override
				public Object resolveEntity(String arg0, String arg1, String arg2,
						String arg3) throws XMLStreamException
					{
					LOG.info("resolveEntity:" +arg0+"/"+arg1+"/"+arg2);
					return null;
					}
				});
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			wx=xof.createXMLEventWriter(System.out, "UTF-8");
			
			if(args.isEmpty())
				{
				LOG.info("Reading from stdin");
				rx=xmlInputFactory.createXMLEventReader(System.in);
				}
			else if(args.size()==1)
				{
				LOG.info("Reading from "+args.get(0));
				rx=xmlInputFactory.createXMLEventReader(IOUtils.openURIForBufferedReading(args.get(0)));
				}
			else
				{
				LOG.error("Illegal number of args");
				return -1;
				}
			run(rx,wx);
			
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(wx);
			CloserUtil.close(rx);
			}

		}
	public static void main(String[] args) {
		new MergeSplittedBlast().instanceMainWithExit(args);
	}
}
