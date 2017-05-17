package com.github.lindenb.jvarkit.tools.blastmapannots;
/**
 * Author:
 * 		Pierre LIndenbaum PhD
 * WWW:
 * 		http://plindenbaum.blogspot.com
 * Mail:
 * 		plindenbaum@yahoo.fr
 * Motivation:
 * 		map the annotations of a genbank uniprot file to a blast hit
 * 		prints the results as a BED file.
 * 
 */
import gov.nih.nlm.ncbi.blast.BlastOutput;
import gov.nih.nlm.ncbi.blast.BlastOutputIterations;
import gov.nih.nlm.ncbi.blast.Hit;
import gov.nih.nlm.ncbi.blast.Hsp;

import gov.nih.nlm.ncbi.blast.Iteration;
import gov.nih.nlm.ncbi.gb.GBFeature;
import gov.nih.nlm.ncbi.gb.GBInterval;
import gov.nih.nlm.ncbi.gb.GBQualifier;
import gov.nih.nlm.ncbi.gb.GBSeq;
import gov.nih.nlm.ncbi.gb.GBSet;


import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Unmarshaller;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.uniprot.Entry;
import org.uniprot.FeatureType;
import org.uniprot.LocationType;
import org.uniprot.Uniprot;
import org.w3c.dom.Document;
import org.xml.sax.EntityResolver;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
/**
BEGIN_DOC

## Example

Download uniprot P04514  ( Rotavirus Non-structural protein 3 )  as <b>XML</b>
```bash
$ curl -o P04514.xml "http://www.uniprot.org/uniprot/P04514.xml"
```
Download the same P04514 as <b>fasta</b>
```bash
$ curl -o P04514.fasta "http://www.uniprot.org/uniprot/P04514.fasta"
```

<b>TblastN</b> P04514.fasta vs a RNA of NSP3 in genbank http://www.ncbi.nlm.nih.gov/nuccore/AY065842.1 and save the ouput as XML:
```xml
<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>tblastn</BlastOutput_program>
(...)
<Hit>
  <Hit_num>1</Hit_num>
  <Hit_id>gi|18139606|gb|AY065842.1|</Hit_id>
  <Hit_def>Rhesus rotavirus nonstructural protein 3 (NSP3) gene, complete cds</Hit_def>
  <Hit_accession>AY065842</Hit_accession>
  <Hit_len>1078</Hit_len>
  <Hit_hsps>
    <Hsp>
      <Hsp_bit-score>546.584</Hsp_bit-score>
      <Hsp_score>1407</Hsp_score>
      <Hsp_evalue>0</Hsp_evalue>
      <Hsp_query-from>1</Hsp_query-from>
      <Hsp_query-to>313</Hsp_query-to>
      <Hsp_hit-from>26</Hsp_hit-from>
      <Hsp_hit-to>964</Hsp_hit-to> <Hsp_qseq>MLKMESTQQMASSIINTSFEAAVVAATSTLELMGIQYDYNEIYTRVKSKFDYVMDDSGVKNNLLGKAATIDQALNGKFGSVMRNKNWMTDSRTVAKLDEDVNKLRMMLSSKGIDQKMRVLNACFSVKRIPGKSSSVIKCTRLMKDKIERGAVEVDDSFVEEKMEVDTVDWKSRYDQLERRFESLKQRVNEKYTTWVQKAKKVNENMYSLQNVISQQQNQIADLQNYCSKLEADLQNKVGSLVSSVEWYLKSMELPDEVKTDIEQQLNSIDTISPINAIDDLEILIRNLIHDYDRTFLMFKGLLRQCNYEYAYE</Hsp_qseq>
      <Hsp_hseq>MLKMESTQQMASSIINSSFEAAVVAATSTLELMGIQYDYNEVYTRVKSKFDLVMDDSGVKNNLIGKAITIDQALNGKFSSAIRNRNWMTDSRTVAKLDEDVNKLRIMLSSKGIDQKMRVLNACFSVKRIPGKSSSIVKCTRLMKDKLERGEVEVDDSFVEEKMEVDTIDWKSRYEQLEKRFESLKHRVNEKYNHWVLKARKVNENMNSLQNVISQQQAHINELQMYNNKLERDLQSKIGSVVSSIEWYLRSMELSDDVKSDIEQQLNSIDQLNPVNAIDDFESILRNLISDYDRLFIMFKGLLQQCNYTYTYE</Hsp_hseq>
      <Hsp_midline>MLKMESTQQMASSIIN SFEAAVVAATSTLELMGIQYDYNE YTRVKSKFD VMDDSGVKNNL GKA TIDQALNGKF S  RN NWMTDSRTVAKLDEDVNKLR MLSSKGIDQKMRVLNACFSVKRIPGKSSS  KCTRLMKDK ERG VEVDDSFVEEKMEVDT DWKSRY QLE RFESLK RVNEKY  WV KA KVNENM SLQNVISQQQ  I  LQ Y  KLE DLQ K GS VSS EWYL SMEL D VK DIEQQLNSID   P NAIDD E   RNLI DYDR F MFKGLL QCNY Y YE</Hsp_midline>
    </Hsp>
  </Hit_hsps>
</Hit>
(...)
</Iteration>
</BlastOutput_iterations>
</BlastOutput>
```
Now produce a BED file with this blast result to map the features of P04514 to AY065842.

```bash
$ java -jar dist/blastmapannots.jar I=P04514.xml B=blast.xml

AY065842	25	961	Non-structural_protein_3	943	+	25961	255,255,255	1	936	25
AY065842	34	469	RNA-binding	970	+	34	469	255,255,255	1	435	34
AY065842	472	640	Dimerization	947	+	472	640	255,255,255	1	168	472
AY065842	532	724	Interaction_with_ZC3H7B	917	+	532	724	255,255,255	1	192	532
AY065842	646	961	Interaction_with_EIF4G1	905	+	646	961	255,255,255	1	315	646
AY065842	520	733	coiled-coil_region	916	+	520	733	255,255,255	1	213	520
```
END_DOC

*/
@Program(name="blastmapannots",
	description="Maps uniprot/genbank annotations on a blast result.",
	keywords={"blast","annotation","genbank","uniprot"}
		)
public class BlastMapAnnotations
	extends Launcher
	{
    private static Logger LOG=Logger.build(BlastMapAnnotations.class).make(); 
	
	private BlastOutput blastOutput=null;
	
    
    @Parameter(names={"-u","-g","--genbank","--uniprot"}, description="XML sequence file Genbank.xml or uniprot.xml.",required=true)
    private File IN=null;
    @Parameter(names={"-APC"}, description="append the sequence accession before the feature name.",required=true)
    private boolean APPEND_ACN=false;	
    @Parameter(names= "--include", description="Restrict to uniprot/feature/type of genbank/feature/key.")
	private Set<String> INCL=new HashSet<String>();
    @Parameter(names= "--exclude", description="Exclude uniprot/feature/type of genbank/feature/key.")
	private Set<String> EXCL=new HashSet<String>();
	

	
	private static class BlastPos
		{
		int query;
		int hit;
		@Override
		public String toString() {
			return "{hit:"+hit+",query:"+query+"}";
			}
		}
	private abstract class Interval
		{
		
		Hit hit;
		Hsp hsp;
		
		protected final BlastPos convertQuery(int qPos1)
			{
			BlastPos left=getHspStart1();
			
			String qS=this.hsp.getHspQseq();
			String hS=this.hsp.getHspHseq();
			for(int i=0;i< qS.length() && i< hS.length();++i)
				{
				if(left.query>=qPos1) break;
				if(isSeqLetter(qS.charAt(i)))
					{
					left.query+=this.queryShift();
					}
				if(isSeqLetter(hS.charAt(i)))
					{
					left.hit+=this.hitShift();
					}
				}
			return left;
			}
		
		protected final int queryShift()
			{
			return 1;
			}
		
		protected final int hitShift()
			{
			int shift=1;
			if(blastOutput.getBlastOutputProgram().equals("tblastn"))
				{
				shift=3;
				}
			else if(blastOutput.getBlastOutputProgram().equals("blastn"))
				{
				shift=1;
				}
			else
				{
				throw new RuntimeException("Sorry program not handled: "+blastOutput.getBlastOutputProgram());
				}
			return shift*(isHspForward()?1:-1);
			}
		
	
		private final boolean isSeqLetter(char c)
			{
			if(c=='-' || c=='*' || Character.isWhitespace(c)) return false;
			if(Character.isLetter(c)) return true;
			throw new IllegalArgumentException("letter: "+c);
			}
		private final boolean isMidMatch(char c)
			{
			if(Character.isWhitespace(c)) return false;
			if(Character.isLetter(c) || c=='|' || c=='+') return true;
			return false;
			}

		
		
		
		protected final BlastPos getHspStart1()
			{
			BlastPos p=new BlastPos();
			p.query=Integer.parseInt(this.hsp.getHspQueryFrom());
			p.hit=Integer.parseInt(this.hsp.getHspHitFrom());
			return p;
			}
		protected final BlastPos getHspEnd1()
			{
			BlastPos p=new BlastPos();
			p.query=Integer.parseInt(this.hsp.getHspQueryTo());
			p.hit=Integer.parseInt(this.hsp.getHspHitTo());
			return p;
			}
		
		protected final boolean isHspForward()
			{
			if(getHspStart1().query>getHspEnd1().query) throw new IllegalStateException();
			return getHspStart1().hit<=getHspEnd1().hit;
			}
		
		public final int getHspQueryStart0()
			{
			return Math.min(getHspStart1().query,getHspEnd1().query)-1;
			}
		public final int getHspQueryEnd0()
			{
			return Math.max(getHspStart1().query,getHspEnd1().query);
			}

		public final int getBedScore()
			{
			int len=featureEnd0()-featureStart0();
			BlastPos left=getHspStart1();
			float match=0f;
			String qS=this.hsp.getHspQseq();
			String hS=this.hsp.getHspHseq();
			String mid=this.hsp.getHspMidline();
			for(int i=0;i< qS.length() && i< hS.length();++i)
				{
				if(left.query-1 >= featureStart0() && left.query-1 < featureEnd0())
					{
					if(isSeqLetter(qS.charAt(i)) )
						{
						match+=(1/3.0);
						}
					if(isSeqLetter(hS.charAt(i)))
						{
						match+=(1/3.0);
						}
					if(isMidMatch(mid.charAt(i)))
						{
						match+=(1/3.0);
						}
					}
				if(isSeqLetter(qS.charAt(i)) )
					{
					left.query+=this.queryShift();
					}
				if(isSeqLetter(hS.charAt(i)))
					{
					left.hit+=this.hitShift();
					}

				}
			int score= (int)((match/len)*1000f);
			if(score<0) score=0;
			if(score>1000)
				{
				LOG.error("SCORE > 10000 in "+toString());
				score=1000;
				}
			if(score==0) LOG.info("score==0 "+match+"/"+len);
			return score;
			}

		public final boolean isFeatureOverlapHsp()
			{
			if(featureEnd0()<=getHspQueryStart0()) return false;
			if(featureStart0()>=getHspQueryEnd0()) return false;
			return true;
			}
		
		public abstract Color getColor();
		
		public final String getBedColor()
			{
			 Color c= getColor();
			 if(c==null) c=Color.WHITE;
			 return ""+c.getRed()+","+c.getGreen()+","+c.getBlue();
			}

		
		
		public final String getChrom()
			{
			return hit.getHitDef();
			}
		protected abstract int featureStart0();
		protected abstract int featureEnd0();

		public abstract String getBedName();
		public abstract int getBedStart();
		public abstract int getBedEnd();
		public abstract char getBedStrand();
		
		private String normalizeName(String s)
			{
			return s.replaceAll("[ \t,]+","_");
			}
		
		public final String toBedString()
			{
			StringBuilder b=new StringBuilder();
			b.append(normalizeName(getChrom())).append('\t');
			b.append(getBedStart()).append('\t');
			b.append(getBedEnd()).append('\t');
			b.append(normalizeName(getBedName())).append('\t');
			b.append(getBedScore()).append('\t');
			b.append(getBedStrand()).append('\t');
			b.append(getBedStart()).append('\t');
			b.append(getBedEnd()).append('\t');
			b.append(getBedColor()).append('\t');
			b.append(1).append('\t');
			b.append(getBedEnd()-getBedStart()).append('\t');
			b.append(0);
			return b.toString();
			}

		
		}
	
	/**===================================================================================
	 * Interval for Uniprot
	 */

	private class UniprotInterval extends Interval
		{
		Entry entry;
		FeatureType featureType;
		
		
		public int getEntryStart1()
			{
			LocationType lt=featureType.getLocation();
			if(lt==null) throw new IllegalStateException();
			if(lt.getPosition()!=null)
				{
				return lt.getPosition().getPosition().intValue();
				}
			else if(lt.getBegin()!=null)
				{
				return lt.getBegin().getPosition().intValue();
				}
			else
				{
				throw new IllegalStateException();
				}
			}
		
		public int getEntryEnd1()
			{
			LocationType lt=featureType.getLocation();
			if(lt==null) throw new IllegalStateException();
			if(lt.getPosition()!=null)
				{
				return lt.getPosition().getPosition().intValue();
				}
			else if(lt.getEnd()!=null)
				{
				return lt.getEnd().getPosition().intValue();
				}
			else
				{
				throw new IllegalStateException();
				}
			}
		
		@Override
		public String toString() {
			return "start1:"+getEntryStart1()+" end1:"+getEntryEnd1()+" acn:"+getBedName()+" start0:"+getEntryStart0()+" end0:"+getEntryEnd0()+"\n"+
					" hsp.start:"+getHspStart1()+" hsp.end:"+getHspEnd1()+" hsp.foward:"+isHspForward()+"\nHSP:overlap-gb:"+isFeatureOverlapHsp();
			}
		
		private String entryName()
			{
			for(String s:entry.getName()) return s;
			for(String s:entry.getAccession()) return s;
			throw new IllegalStateException();
			}
		
		@Override 
		public String getBedName()
			{
			String s=this.featureType.getType();
			if(this.featureType.getDescription()!=null && !s.equals("sequence conflict"))
				{
				s= this.featureType.getDescription();
				}
			if(s.equals("sequence conflict") && featureType.getOriginal()!=null && featureType.getVariation()!=null && !featureType.getVariation().isEmpty())
				{
				s+=":"+featureType.getOriginal()+"/";
				for(String v:featureType.getVariation())
					{
					if(!s.endsWith("/")) s+=",";
					s+=v;
					}
				}
			if(!APPEND_ACN) return s;
			return entryName()+":"+s;
			}
		
		private int getEntryStart0()
			{
			return getEntryStart1()-1;
			}
		private int getEntryEnd0()
			{
			return getEntryEnd1();
			}
		
		@Override
		protected int featureEnd0()
			{
			return getEntryEnd1();
			}
		@Override
		protected int featureStart0()
			{
			return getEntryStart0();
			}
		
		
		@Override
		public int getBedStart()
			{
			return Math.min(convertQuery(getEntryStart1()).hit,convertQuery(getEntryEnd1()).hit)-1;
			}
		
		@Override
		public int getBedEnd()
			{
			int bedEnd= Math.max(convertQuery(getEntryStart1()).hit,convertQuery(getEntryEnd1()).hit);
			if(bedEnd==getBedStart())
				{
				LOG.warning("Empty bed:");
				}
			return bedEnd;
			}

		
		@Override
		public char getBedStrand()
			{
			int str=1;
			if(!isHspForward()) str*=-1;
			return str==1?'+':'-';
			}
		
		@Override
		public Color getColor()
			{
			String fkey="";//TODO
			if(fkey.equals("CDS"))
				{
				return Color.YELLOW;
				}
			if(fkey.equals("gene"))
				{
				return Color.ORANGE;
				}
			if(fkey.equals("gene"))
				{
				return Color.GREEN;
				}
			return Color.WHITE;
			}
		
		}
	static enum QualKey
		{
		region_name,product,mol_type,gene,locus_tag,site_type,note
		};
	/**===================================================================================
	 * Interval for Genbank
	 */
	private class GenbankInterval extends Interval
		{
		
		GBSeq gbSeq;
		GBInterval gbInterval;
		GBFeature gbFeature;

		

		
		public int getGBStart1()
			{
			if(gbInterval.getGBIntervalPoint()!=null)
				{
				return Integer.parseInt(gbInterval.getGBIntervalPoint());
				}
			else if(gbInterval.getGBIntervalFrom()!=null)
				{
				return Integer.parseInt(gbInterval.getGBIntervalFrom());
				}
			else
				{
				throw new IllegalStateException();
				}
			}
		
		public int getGBEnd1()
			{
			if(gbInterval.getGBIntervalPoint()!=null)
				{
				return Integer.parseInt(gbInterval.getGBIntervalPoint());
				}
			else if(gbInterval.getGBIntervalTo()!=null)
				{
				return Integer.parseInt(gbInterval.getGBIntervalTo());
				}
			else
				{
				throw new IllegalStateException();
				}
			}
		
		@Override
		public String toString() {
			return "start1:"+getGBStart1()+" end1:"+getGBEnd1()+" forward:"+isGbForward()+" acn:"+getBedName()+" start0:"+getGBStart0()+" end0:"+getGBEnd0()+"\n"+
					" hsp.start:"+getHspStart1()+" hsp.end:"+getHspEnd1()+" hsp.foward:"+isHspForward()+"\nHSP:overlap-gb:"+isFeatureOverlapHsp();
			}
		
		
		private String gbName()
			{
			String s=null;
			s=gbSeq.getGBSeqLocus();
			if(s!=null) return s;
			s=gbSeq.getGBSeqAccessionVersion();
			if(s!=null) return s;
			s=gbSeq.getGBSeqEntryVersion();
			if(s!=null) return s;
			return "?";
			}
		
		@Override
		public String getBedName()
			{
			String fkey=this.gbFeature.getGBFeatureKey();
			TreeMap<QualKey, String> t=new TreeMap<QualKey, String>();
			
			if(this.gbFeature.getGBFeatureQuals()!=null &&
					this.gbFeature.getGBFeatureQuals().getGBQualifier()!=null
					)
				{
				for(GBQualifier q:this.gbFeature.getGBFeatureQuals().getGBQualifier())
					{
					String key=q.getGBQualifierName();
					for(QualKey qk:QualKey.values())
						{
						if(qk.name().equals(key))
							{
							t.put(qk, q.getGBQualifierValue());
							break;
							}
						}
					
					}
				}
			if(t.isEmpty())
				{
				LOG.info("not qual for "+fkey);
				if(!APPEND_ACN) return fkey;
				return gbName()+":"+fkey;
				}
			if(!APPEND_ACN) return fkey+":"+t.get(t.keySet().iterator().next());
			return gbName()+":"+fkey+":"+t.get(t.keySet().iterator().next());
			}
		
		public boolean isGbForward()
			{
			return getGBStart1()<=getGBEnd1();
			}
		
		private int getGBStart0()
			{
			return Math.min(getGBStart1(), getGBEnd1())-1;
			}
		private int getGBEnd0()
			{
			return Math.max(getGBStart1(), getGBEnd1());
			}
		
		
		@Override
		protected int featureEnd0()
			{
			return getGBEnd0();
			}
		@Override
		protected int featureStart0()
			{
			return getGBStart0();
			}
		
	
		
		
		

		@Override
		public int getBedStart()
			{
			return Math.min(convertQuery(getGBStart1()).hit,convertQuery(getGBEnd1()).hit)-1;
			}
		@Override
		public int getBedEnd()
			{
			return Math.max(convertQuery(getGBStart1()).hit,convertQuery(getGBEnd1()).hit);
			}

		
		
		@Override
		public char getBedStrand()
			{
			int str=1;
			if(!isGbForward()) str*=-1;
			if(!isHspForward()) str*=-1;
			return str==1?'+':'-';
			}
		@Override
		public Color getColor()
			{
			String fkey=this.gbFeature.getGBFeatureKey();
			if(fkey.equals("CDS"))
				{
				return Color.YELLOW;
				}
			if(fkey.equals("gene"))
				{
				return Color.ORANGE;
				}
			if(fkey.equals("gene"))
				{
				return Color.GREEN;
				}
			return Color.WHITE;
			}
		}

	private void printUniprot(Uniprot uniprotSet)
		{
		
		if(uniprotSet.getEntry().isEmpty())
			{
			LOG.warn("empty uniprot entry.");
			return;
			}
		
		if(uniprotSet.getEntry().size()>1)
			{
			LOG.warn("entry contains more than one sequence.");
			}
		for(Entry entry:uniprotSet.getEntry())
			{	
			BlastOutputIterations iterations=this.blastOutput.getBlastOutputIterations();
			for(Iteration iteration:iterations.getIteration())
				{
				for(FeatureType feature:entry.getFeature())
					{
					if(!acceptfeature(feature.getType())) continue;

					
					for(Hit hit:iteration.getIterationHits().getHit())
						{
						for(Hsp hsp :hit.getHitHsps().getHsp())
							{
							UniprotInterval bi=new UniprotInterval();
						
							bi.entry=entry;
							bi.featureType=feature;
							bi.hit=hit;
							bi.hsp=hsp;
							LOG.debug("interval "+bi);
							if(!bi.isFeatureOverlapHsp())
								{
								continue;
								}
							stdout().println(bi.toBedString());
							}
						}
					
					}
				break;
				}
			break;
			}
		
		//System.err.println("OK");
		
		}	
	
	
	
	private void printGB(GBSet gbSet)
		{
		for(GBSeq gbSeq:gbSet.getGBSeq())
			{	
			BlastOutputIterations iterations=this.blastOutput.getBlastOutputIterations();
			for(Iteration iteration:iterations.getIteration())
				{
				for(GBFeature feature:gbSeq.getGBSeqFeatureTable().getGBFeature())
					{
					if(feature.getGBFeatureIntervals()==null) continue;
					
					if(!acceptfeature(feature.getGBFeatureKey())) continue;
					
					

					
					for(GBInterval interval:feature.getGBFeatureIntervals().getGBInterval())
						{
		
						for(Hit hit:iteration.getIterationHits().getHit())
							{
							for(Hsp hsp :hit.getHitHsps().getHsp())
								{
								GenbankInterval bi=new GenbankInterval();
								bi.gbSeq=gbSeq;
								bi.gbFeature=feature;
								bi.gbInterval=interval;
								bi.hit=hit;
								bi.hsp=hsp;
								LOG.debug("interval "+bi);
								if(!bi.isGbForward()) LOG.info("CHECK INTERVAL REVERSE");
								if(!bi.isFeatureOverlapHsp()) continue;
								stdout().println(bi.toBedString());
								}
							}
						
						}
					}
				break;
				}
			}
		
		//System.err.println("OK");
		
		}
	
	@Override
	public int doWork(List<String> args) {
		try
			{
			/** xml parser */
			 DocumentBuilder docBuilder;
			/** transforms XML/DOM to GBC entry */
			 Unmarshaller unmarshaller;
	
			//create a DOM parser
			DocumentBuilderFactory f=DocumentBuilderFactory.newInstance();
			f.setCoalescing(true);
			//f.setNamespaceAware(true); no, why does it break the parsing of uniprot ??
			f.setValidating(false);
			f.setExpandEntityReferences(true);
			docBuilder= f.newDocumentBuilder();
			docBuilder.setEntityResolver(new EntityResolver()
				{
				@Override
				public InputSource resolveEntity(String publicId, String systemId)
						throws SAXException, IOException
					{
					return new InputSource(new StringReader(""));
					}
				});
			//create a Unmarshaller for NCBI
			JAXBContext jc = JAXBContext.newInstance("gov.nih.nlm.ncbi.gb:gov.nih.nlm.ncbi.blast:org.uniprot");
			unmarshaller=jc.createUnmarshaller();
	
			
			
			LOG.info("reading entry "+IN);
			Document domEntry=docBuilder.parse(IN);
			GBSet gbSet=null;
			Uniprot uniprotSet=null;
			if("GBSet".equals(domEntry.getDocumentElement().getNodeName()))
				{
				LOG.info("parsing as GBSet");
				gbSet=unmarshaller.unmarshal(domEntry,GBSet.class).getValue();	
				}
			else if("uniprot".equals(domEntry.getDocumentElement().getNodeName()))
				{
				LOG.info("parsing as Uniprot "+domEntry.getDocumentElement());
				uniprotSet=unmarshaller.unmarshal(domEntry,Uniprot.class).getValue();	
				//LOG.info(uniprotSet.getEntry().size());
				//jc.createMarshaller().marshal(uniprotSet, System.err); 
				}
			else
				{
				LOG.info("unknown root element:"+domEntry.getDocumentElement().getNodeName());
				return -1;
				}
			Document blastDom;
			if(args.size()==1)
				{
				LOG.info("reading "+args.get(0));
				blastDom=docBuilder.parse(new File(args.get(0)));
				}
			else if(args.isEmpty())
				{
				LOG.info("reading from stdin");
				blastDom=docBuilder.parse(stdin());
				}
			else
				{
				LOG.error("Illegal number of args");
				return -1;
				}
			this.blastOutput=unmarshaller.unmarshal(blastDom,BlastOutput.class).getValue();	
			if(uniprotSet!=null) printUniprot(uniprotSet);
			if(gbSet!=null) printGB(gbSet);
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}	
		}
	private boolean acceptfeature(String s)
		{
		if(s==null || s.isEmpty()) return false;
		if(!this.INCL.isEmpty())
			{
			if(!INCL.contains(s))
				{
				return false;
				}
			}
		if(!this.EXCL.isEmpty())
			{
			if(EXCL.contains(s))
				{
				return false;
				}
			}
		return true;
		}
	
	
	public static void main(String[] args)
		{
		/* force javac to compile those */
		@SuppressWarnings("unused")
		gov.nih.nlm.ncbi.blast.ObjectFactory of1=null;
		@SuppressWarnings("unused")
		gov.nih.nlm.ncbi.gb.ObjectFactory of2=null;
		@SuppressWarnings("unused")
		org.uniprot.ObjectFactory of3=null;
		new BlastMapAnnotations().instanceMainWithExit(args);
		}
	}
