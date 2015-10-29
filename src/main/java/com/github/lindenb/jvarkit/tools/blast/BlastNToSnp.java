/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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

*/
package com.github.lindenb.jvarkit.tools.blast;

import java.io.PrintWriter;
import java.util.Collection;

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

import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.command.Command;



public class BlastNToSnp extends AbstractBlastNToSnp
	{
	/* force javac to compile */
	@SuppressWarnings("unused")
	private gov.nih.nlm.ncbi.blast.ObjectFactory _ignore_for_javac=null;

	
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(BlastNToSnp.class);

	
	@Override
	public Command createCommand() {
		return new MyCommand();
		}

	static private class MyCommand extends AbstractBlastNToSnp.AbstractBlastNToSnpCommand
		{    
		private Unmarshaller unmarshaller;
		private Marshaller marshaller;

	
	
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
		PrintWriter pw,
		XMLEventReader r
		)
		throws XMLStreamException,JAXBException
	{
	 long numIterations=0L;
	 long numPerfectMath=0L;
	 long numNoHit=0L;
	 
	
	pw.print("#query");
	pw.print('\t');
	pw.print("hit");
	pw.print('\t');
	pw.print(("hit-index"));
	pw.print('\t');
	pw.print("hsp-index");
	
	pw.print('\t');
	pw.print("query-POS");
	pw.print('\t');
	pw.print("hit-POS");
	pw.print('\t');
	pw.print("STRAND");
	
	pw.print('\t');
	pw.print("REF(hit)");
	pw.print('\t');
	pw.print("ALT(query)");
	pw.print('\t');
	pw.print("blast.align_length");
	pw.print('\t');
	pw.print("blast.hit.var");
	pw.print('\t');
	pw.print("blast.query.var");
	pw.print('\t');
	pw.print("blast.mid.var");
	
	pw.println();
	for(;;)
		{

		Iteration iter1=peekIteration(r);
		if(iter1==null) break;
		++numIterations;
		if(iter1.getIterationHits().getHit().isEmpty())
			{
			LOG.info("No hit found for "+iter1.getIterationQueryDef());
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
						if(k==-1 || k-j> minGapSize) break;
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
						
					pw.print(iter1.getIterationQueryDef());
					pw.print('\t');
					pw.print(hit.getHitDef());
					pw.print('\t');
					pw.print((1+hit_loop));
					pw.print('\t');
					pw.print((1+hsp_loop));
					
					pw.print('\t');
					pw.print(query_index);
					pw.print('\t');
					pw.print(hit_index);
					pw.print('\t');
					pw.print(hit_shift==1?'+':'-');
					
					pw.print('\t');
					pw.print(ref);
					pw.print('\t');
					pw.print(alt);
					
					
					pw.print('\t');
					pw.print(align_length);
					pw.print('\t');
					pw.print(hspHSeq.substring(i, j).replace(' ', '.'));
					pw.print('\t');
					pw.print(hspQseq.substring(i, j).replace(' ', '.'));
					pw.print('\t');
					pw.print(hspMid.substring(i, j).replace(' ', '.'));

					
					pw.println();
					
					//marshaller.marshal(new JAXBElement<Hsp>(new QName("Hsp"), Hsp.class, hsp), System.out);
					//pw.println();
					
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
				if(pw.checkError()) break;
				}
			
			}//end loop Hit
		
		if(!found_mismatch)
			{
			numPerfectMath++;
			}
		
		}//end while read
	LOG.info("ITERATIONS : " +numIterations);
	LOG.info("ONLY_PERFECT_MATCH : " +numPerfectMath);
	LOG.info("NO_HIT : " +numNoHit);

	}

	@Override
	protected Collection<Throwable> call(String inputName) throws Exception
			{
			PrintWriter pw=null;
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
						LOG.info("resolveEntity:" +arg0+"/"+arg1+"/"+arg2);
						return null;
						}
					});
				if(inputName==null)
					{
					LOG.info("Reading from stdin");
					rx=xmlInputFactory.createXMLEventReader(stdin());
					}
				else
					{
					LOG.info("Reading from "+inputName);
					rx=xmlInputFactory.createXMLEventReader(IOUtils.openURIForBufferedReading(inputName));
					}
				pw = openFileOrStdoutAsPrintWriter();
				run(pw,rx);
				pw.flush();
				return RETURN_OK;
				}
			catch(Exception err)
				{
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(rx);
				CloserUtil.close(pw);
				}
			}
		}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new BlastNToSnp().instanceMainWithExit(args);

	}

}
