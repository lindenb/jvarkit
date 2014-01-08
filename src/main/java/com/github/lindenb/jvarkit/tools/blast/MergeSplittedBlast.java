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

import net.sf.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.bio.blast.BlastHspAlignment;

public class MergeSplittedBlast extends AbstractCommandLineProgram {
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
			error(err,"expect (chrom:start-end:length) Cannot parse Hit.def: "+s );
			return null;
			}
		return split;
		}
	
	
	private boolean overlap(int a_start1,int a_end1,int b_start1,int b_end1)
		{
		info(" overlapp ["+a_start1+"-"+a_end1+"] vs ["+b_start1+"-"+b_end1+"]");
		if(b_end1+1 < a_start1) return false;
		if(a_end1+1 < b_start1) return false;
		info(" ok overlap");
		return true;
		}
	private Hsp merge(Hsp hh0,Hsp hh1)
		{
		info("ICIIII");
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
			debug("strand");
			return null;
			}
		
		
		if(!overlap(aln0.getHitFrom1(),aln0.getHitTo1(),aln1.getHitFrom1(),aln1.getHitTo1()))
			{
			
			return null;
			}
		if(!overlap(aln0.getQueryFrom1(),aln0.getQueryTo1(),aln1.getQueryFrom1(),aln1.getQueryTo1()))
			{
			return null;
			}
		//hit1 is contained in hit0
		if(aln0.getQueryFrom1() <= aln1.getQueryFrom1() && aln1.getQueryTo1()<=aln0.getQueryTo1())
			{
			return hh0;
			}
		
		StringBuilder qsb=new StringBuilder();
		StringBuilder msb=new StringBuilder();
		StringBuilder hsb=new StringBuilder();
		for(BlastHspAlignment.Align a:aln0)
			{
			if(a.getQueryIndex1()>=aln1.getQueryFrom1()) break;
			qsb.append(a.getQueryChar());
			msb.append(a.getMidChar());
			hsb.append(a.getHitChar());
			}
		for(BlastHspAlignment.Align a:aln1)
			{
			if(a.getQueryIndex1()<aln1.getQueryFrom1()) continue;
			qsb.append(a.getQueryChar());
			msb.append(a.getMidChar());
			hsb.append(a.getHitChar());
			}
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
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/MergeSplittedBlast";
		}
	@Override
	public String getProgramDescription() {
		return "merge blast Hits from splitted BLAST database. See http://www.biostars.org/p/90186/";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+""))!=-1)
			{
			switch(c)
				{
				default:
					{
					switch(handleOtherOptions(c, opt))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		
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
					info("resolveEntity:" +arg0+"/"+arg1+"/"+arg2);
					return null;
					}
				});
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			wx=xof.createXMLEventWriter(System.out, "UTF-8");
			
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				rx=xmlInputFactory.createXMLEventReader(System.in);
				}
			else if(opt.getOptInd()+1==args.length)
				{
				info("Reading from "+args[opt.getOptInd()]);
				rx=xmlInputFactory.createXMLEventReader(IOUtils.openURIForBufferedReading(args[opt.getOptInd()]));
				}
			else
				{
				error("Illegal number of args");
				return -1;
				}
			run(rx,wx);
			
			return 0;
			}
		catch(Exception err)
			{
			error(err);
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
