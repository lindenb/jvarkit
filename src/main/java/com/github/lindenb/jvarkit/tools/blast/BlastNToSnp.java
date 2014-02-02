package com.github.lindenb.jvarkit.tools.blast;

import gov.nih.nlm.ncbi.blast.Hit;
import gov.nih.nlm.ncbi.blast.Hsp;
import gov.nih.nlm.ncbi.blast.Iteration;


import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLResolver;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.XMLEvent;
import net.sf.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;



public class BlastNToSnp extends AbstractCommandLineProgram {
	private Unmarshaller unmarshaller;
	private Marshaller marshaller;
	
	/* force javac to compile */
	@SuppressWarnings("unused")
	private gov.nih.nlm.ncbi.blast.ObjectFactory _ignore_for_javac=null;
	
	private BlastNToSnp()
		{
		
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/BlastnToSnp";
		}
	
	@Override
	public String getProgramDescription() {
		return "print indel/mismatch in a blastn stream";
		}
	
	private Iteration peekIteration(XMLEventReader r) throws XMLStreamException,JAXBException
		{
		while(r.hasNext())
			{
			XMLEvent evt=r.peek();
			if(!(evt.isStartElement()))
				{
				r.next();
				continue;
				}
			String name=evt.asStartElement().getName().getLocalPart();
			if(name.equals("BlastOutput_program"))
				{
				r.next();
				String prog=r.getElementText();
				if(!"blastn".equals(prog))
					{
					throw new XMLStreamException("Not a blastn input:"+prog);
					}
				continue;
				}
			
			if(!name.equals("Iteration"))
				{
				r.next();
				continue;
				}
			
			return this.unmarshaller.unmarshal(r, Iteration.class).getValue();
			}
		return null;
		}
	
	
private void run(
		XMLEventReader r
		)
		throws XMLStreamException,JAXBException
	{
	 long numIterations=0L;
	 long numPerfectMath=0L;
	 long numNoHit=0L;

	
	System.out.print("#query");
	System.out.print('\t');
	System.out.print("hit");
	System.out.print('\t');
	System.out.print(("hit-index"));
	System.out.print('\t');
	System.out.print("hsp-index");
	
	System.out.print('\t');
	System.out.print("query-POS");
	System.out.print('\t');
	System.out.print("hit-POS");
	System.out.print('\t');
	System.out.print("STRAND");
	
	System.out.print('\t');
	System.out.print("REF(hit)");
	System.out.print('\t');
	System.out.print("ALT(query)");
	System.out.print('\t');
	System.out.print("blast.align_length");
	System.out.print('\t');
	System.out.print("blast.hit.var");
	System.out.print('\t');
	System.out.print("blast.query.var");
	System.out.print('\t');
	System.out.print("blast.mid.var");
	
	System.out.println();
	for(;;)
		{

		Iteration iter1=peekIteration(r);
		if(iter1==null) break;
		++numIterations;
		if(iter1.getIterationHits().getHit().isEmpty())
			{
			info("No hit found for "+iter1.getIterationQueryDef());
			++numNoHit;
			continue;
			}
		boolean found_mismatch=false;
		for(int hit_loop=0; hit_loop <  iter1.getIterationHits().getHit().size();++hit_loop)
			{
			final Hit hit=iter1.getIterationHits().getHit().get(hit_loop);
			for(int  hsp_loop =0;hsp_loop< hit.getHitHsps().getHsp().size();++hsp_loop)
				{
				final Hsp hsp = hit.getHitHsps().getHsp().get(hsp_loop);
				String hspQseq=hsp.getHspQseq();
				String hspMid=hsp.getHspMidline();
				String hspHSeq=hsp.getHspHseq();
				
				int hsp_query_from=Integer.parseInt(hsp.getHspQueryFrom());
				int hsp_query_to=Integer.parseInt(hsp.getHspQueryTo());
				int hsp_hit_from=Integer.parseInt(hsp.getHspHitFrom());
				int hsp_hit_to=Integer.parseInt(hsp.getHspHitTo());
				final int align_length=Integer.parseInt(hsp.getHspAlignLen());
				final int hit_shift=(hsp_hit_from>hsp_hit_to?-1:1);
				
				int i=0;
				int query_index = hsp_query_from;
				int hit_index = hsp_hit_from;
				
				while( i< align_length )
					{
					char ch= hspHSeq.charAt(i);
					char cq= hspQseq.charAt(i);
					char cm= hspMid.charAt(i);
					if(cm=='|')
						{
						++i;
						query_index++;
						hit_index+=hit_shift;
						continue;
						}
					found_mismatch=true;
					
					int j=i+1;
					for(;;)
						{
						int k=hspMid.indexOf(' ', j);
						if(k==-1 || k-j>3) break;
						j=k+1;
						}
					
					
					String ref=new String(hspHSeq.substring(i, j)).replaceAll("[\\- ]", "");
					String alt=new String(hspQseq.substring(i, j)).replaceAll("[\\- ]", "");
					
					
					
					if(hit_shift<0)
						{
						StringBuilder sb=new StringBuilder(alt.length());
						for(int x=alt.length()-1;x>=0;--x)
							{
							sb.append(AcidNucleics.complement(alt.charAt(x)));
							}
						alt=sb.toString();
						
						sb=new StringBuilder(ref.length());
						for(int x=ref.length()-1;x>=0;--x)
							{
							sb.append(AcidNucleics.complement(ref.charAt(x)));
							}
						ref=sb.toString();
						
						
						}
						
					System.out.print(iter1.getIterationQueryDef());
					System.out.print('\t');
					System.out.print(hit.getHitDef());
					System.out.print('\t');
					System.out.print((1+hit_loop));
					System.out.print('\t');
					System.out.print((1+hsp_loop));
					
					System.out.print('\t');
					System.out.print(query_index);
					System.out.print('\t');
					System.out.print(hit_index);
					System.out.print('\t');
					System.out.print(hit_shift==1?'+':'-');
					
					System.out.print('\t');
					System.out.print(ref);
					System.out.print('\t');
					System.out.print(alt);
					
					
					System.out.print('\t');
					System.out.print(align_length);
					System.out.print('\t');
					System.out.print(hspHSeq.substring(i, j).replace(' ', '.'));
					System.out.print('\t');
					System.out.print(hspQseq.substring(i, j).replace(' ', '.'));
					System.out.print('\t');
					System.out.print(hspMid.substring(i, j).replace(' ', '.'));

					
					System.out.println();
					
					//marshaller.marshal(new JAXBElement<Hsp>(new QName("Hsp"), Hsp.class, hsp), System.out);
					//System.out.println();
					
					while(i<j)
						{
						ch= hspHSeq.charAt(i);
						cq= hspQseq.charAt(i);

						if(ch!='-' && ch!=' ')
							{
							hit_index+=hit_shift;
							}
						if(cq!='-' && cq!=' ')
							{
							query_index++;
							}
						++i;
						}
					}
				if(hit_index-hit_shift!=hsp_hit_to)
					{
					marshaller.marshal(new JAXBElement<Hsp>(new QName("Hsp"), Hsp.class, hsp), System.err);
					throw new IllegalStateException("boum "+hit_index+" vs "+hsp_hit_to);
					}
				if(query_index-1!=hsp_query_to)
					{
					throw new IllegalStateException("boum "+query_index+" vs "+hsp_query_to);
					}
				if(System.out.checkError()) break;
				}
			
			}//end loop Hit
		
		if(!found_mismatch)
			{
			numPerfectMath++;
			}
		
		}//end while read
	info("ITERATIONS : " +numIterations);
	info("ONLY_PERFECT_MATCH : " +numPerfectMath);
	info("NO_HIT : " +numNoHit);

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
		while((c=opt.getopt(args,getGetOptDefault()))!=-1)
			{
			switch(c)
				{
			
				default:
					{
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		
		XMLEventReader rx=null;
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
			run(rx);
			
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(rx);
			}
		}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new BlastNToSnp().instanceMainWithExit(args);

	}

}
