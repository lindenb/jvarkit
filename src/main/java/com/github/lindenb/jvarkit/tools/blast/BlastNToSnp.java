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

*/
package com.github.lindenb.jvarkit.tools.blast;

import java.io.File;
import java.io.PrintWriter;
import java.util.List;

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

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**

BEGIN_DOC

### Example

```
$  java -jar dist/blastn2snp.jar < blast.xml | column -t
```

```
#query              hit                                                                                              hit-index  hsp-index  query-POS  hit-POS    STRAND  REF(hit)  ALT(query)  blast.align_length  blast.hit.var  blast.query.var  blast.mid.var
No definition line  Homo sapiens chromosome 6, alternate assembly CHM1_1.1                                           1          9          21         74567818   -       T         A           18                  A              T                .
No definition line  Homo sapiens chromosome 6, alternate assembly HuRef                                              2          9          21         71600901   -       T         A           18                  A              T                .
No definition line  Homo sapiens chromosome 6, GRCh37.p13 Primary Assembly                                           3          9          21         74401398   -       T         A           18                  A              T                .
No definition line  Homo sapiens chromosome 5, alternate assembly CHM1_1.1                                           4          1          7          107821121  -       A         G           28                  T              C                .
No definition line  Homo sapiens chromosome 5, alternate assembly CHM1_1.1                                           4          9          16         14262358   +       G         C           18                  G              C                .
No definition line  Homo sapiens chromosome 5, alternate assembly CHM1_1.1                                           4          13         8          132662461  -       T         C           18                  A              G                .
No definition line  Homo sapiens chromosome 5, alternate assembly CHM1_1.1                                           4          20         14         170329095  -       G         C           18                  C              G                .
No definition line  Homo sapiens chromosome 5, alternate assembly HuRef                                              5          1          7          103561224  -       A         G           28                  T              C                .
No definition line  Homo sapiens chromosome 5, alternate assembly HuRef                                              5          9          16         14234054   +       G         C           18                  G              C                .
No definition line  Homo sapiens chromosome 5, alternate assembly HuRef                                              5          13         8          128416747  -       T         C           18                  A              G                .
No definition line  Homo sapiens chromosome 5, alternate assembly HuRef                                              5          19         14         165993804  -       G         C           18                  C              G                .
No definition line  Homo sapiens chromosome 5, GRCh37.p13 Primary Assembly                                           6          1          7          108388040  -       A         G           28                  T              C                .
No definition line  Homo sapiens chromosome 5, GRCh37.p13 Primary Assembly                                           6          9          16         14262514   +       G         C           18                  G              C                .
No definition line  Homo sapiens chromosome 5, GRCh37.p13 Primary Assembly                                           6          13         8          133231874  -       T         C           18                  A              G                .
No definition line  Homo sapiens chromosome 5, GRCh37.p13 Primary Assembly                                           6          19         14         170896537  -       G         C           18                  C              G                .
No definition line  Homo sapiens chromosome 19, alternate assembly CHM1_1.1                                          8          13         18         54878835   -       A         C           18                  T              G                .
No definition line  Homo sapiens chromosome 19, alternate assembly CHM1_1.1                                          8          14         18         54999463   +       T         G           18                  T              G                .
No definition line  Homo sapiens chromosome 19, alternate assembly CHM1_1.1                                          8          15         21         55131801   -       G         A           18                  C              T                .
No definition line  Homo sapiens chromosome 19, alternate assembly HuRef                                             9          14         18         51207279   -       A         C           18                  T              G                .
No definition line  Homo sapiens chromosome 19, alternate assembly HuRef                                             9          15         18         51329261   +       T         G           18                  T              G                .
No definition line  Homo sapiens chromosome 19, alternate assembly HuRef                                             9          16         21         51461318   -       G         A           18                  C              T                .

```

END_DOC
 */
@Program(name="blastn2snp",
	keywords={"blast","snp"},
	description="print indel/mismatch in a blastn stream",
	biostars={89151}
	)
public class BlastNToSnp extends Launcher
{
	private static final Logger LOG = Logger.build(BlastNToSnp.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-n","--gapsize"},description="min gap")
	private int minGapSize = 3 ;

	private Unmarshaller unmarshaller;
	private Marshaller marshaller;
	/* force javac to compile */
	@SuppressWarnings("unused")
	private gov.nih.nlm.ncbi.blast.ObjectFactory _ignore_for_javac=null;
	
	
	private Iteration peekIteration(final XMLEventReader r) throws XMLStreamException,JAXBException
		{
		while(r.hasNext())
			{
			final XMLEvent evt=r.peek();
			if(!(evt.isStartElement()))
				{
				r.next();
				continue;
				}
			final String name=evt.asStartElement().getName().getLocalPart();
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
					marshaller.marshal(new JAXBElement<Hsp>(new QName("Hsp"), Hsp.class, hsp),stderr());
					throw new IllegalStateException("Error expected hit-index= "+hit_index+"-"+hit_shift+" == hsp-hit-to="+hsp_hit_to);
					}
				if(query_index-1!=hsp_query_to)
					{
					marshaller.marshal(new JAXBElement<Hit>(new QName("Hit"), Hit.class, hit),stderr());
					throw new IllegalStateException("query_index "+query_index+"(1 != hsp_query_to:"+hsp_query_to);
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
	public int doWork(final List<String> args) {
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
			final String inputName= oneFileOrNull(args);
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
			pw = super.openFileOrStdoutAsPrintWriter(outputFile);
			run(pw,rx);
			pw.flush();
			return RETURN_OK;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(rx);
			CloserUtil.close(pw);
			}
		}

	/**
	 * @param args
	 */
	public static void main(final String[] args) {
		new BlastNToSnp().instanceMainWithExit(args);

	}

}
