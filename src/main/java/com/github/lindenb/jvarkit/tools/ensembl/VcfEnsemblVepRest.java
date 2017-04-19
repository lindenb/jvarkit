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
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.ensembl;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Unmarshaller;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import javax.xml.transform.stream.StreamSource;

import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpPost;
import org.apache.http.entity.ByteArrayEntity;
import org.apache.http.entity.ContentType;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;
import org.ensembl.vep.*;
import org.w3c.dom.Attr;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.TeeInputStream;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

@Program(name="vcfensemblvep",description="Annotate a VCF with ensembl REST API")
public class VcfEnsemblVepRest 
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfEnsemblVepRest.class).make();

	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;


	@Parameter(names={"-s","--server"},description="REST server")
	private String server = "http://grch37.rest.ensembl.org";

	@Parameter(names={"-e","--extension"},description="Path extension")
	private String extension = "/vep/homo_sapiens/region";

	@Parameter(names={"-n","--batchSize"},description="batch size")
	private int batchSize = 100 ;

	@Parameter(names={"-x","--base64"},description="save whole XML document as xml base 64")
	private boolean xmlBase64 = false;

	@Parameter(names={"-T","--tee"},description="'Tee' xml response to stderr")
	private boolean teeResponse = false;

	public static final String TAG="VEPTRCSQ";
	@SuppressWarnings("unused")
	private static final ObjectFactory _fool_javac=null;
	private Unmarshaller unmarshaller=null;
	private DocumentBuilder documentBuilder;
	private Transformer xmlSerializer;
	private CloseableHttpClient httpClient = null;
	
	
	private static String createInputContext(VariantContext ctx)
		{
		StringBuilder sb=new StringBuilder();
		sb.append(ctx.getContig()).
			append(" ").
			append(ctx.getStart()).
			append(" ").
			append(!ctx.hasID()?".":ctx.getID()).
			append(" ").
			append(ctx.getReference().getBaseString()).
			append(" ")
			;
		List<Allele> alts=ctx.getAlternateAlleles();
		if(alts.isEmpty())
			{
			sb.append(".");
			}
		else
			{
			for(int j=0;j< alts.size();++j )
				{
				if(j>0) sb.append(",");
				sb.append(alts.get(j).getBaseString());
				}
			}
		sb.append(" . . .");
		return sb.toString();
		}
	
	private static String empty(Object s)
		{
		return s==null || String.valueOf(s).trim().isEmpty()?"":String.valueOf(s);
		}
	
	private long lastMillisec=-1L;
	private Object generic_vep(List<VariantContext> contexts,boolean xml_answer) throws IOException
		{
		LOG.info("Running VEP "+contexts.size());
		InputStream response =null;
		javax.xml.transform.Source inputSource=null;
		HttpPost httpPost = null;
		try {
		    if ( this.lastMillisec!=-1L && this.lastMillisec+ 5000<  System.currentTimeMillis())
		    	{
		    	try {Thread.sleep(1000);} catch(Exception err){}
		    	}
				 
		    httpPost = new HttpPost(this.server + this.extension);
			 
			 
			 StringBuilder queryb=new StringBuilder();
			 queryb.append("{ \"variants\" : [");
			 for(int i=0;i< contexts.size();++i)
			 	{
				VariantContext ctx=contexts.get(i);
				if(i>0) queryb.append(",");
				queryb.append("\"").
					append(createInputContext(ctx)).
					append("\"");
			 	}
			 queryb.append("]");
			 for(String s: new String[]{"canonical","ccds","domains","hgvs","numbers","protein","xref_refseq"})
			 	{
				 queryb.append(",\"").append(s).append("\":1");
			 	}
			 queryb.append("}");
			 byte postBody[] = queryb.toString().getBytes();

			 httpPost.setHeader("Content-Type",ContentType.APPLICATION_JSON.getMimeType());
			 httpPost.setHeader("Accept",ContentType.TEXT_XML.getMimeType());
			 //httpPost.setHeader("Content-Length", Integer.toString(postBody.length));
			 httpPost.setEntity(new ByteArrayEntity(postBody, ContentType.APPLICATION_JSON));
			 
			 
			 
			 final CloseableHttpResponse httpResponse = httpClient.execute(httpPost);
			 
			 int responseCode = httpResponse.getStatusLine().getStatusCode();
			  
			 if(responseCode != 200)
			 	{
				throw new RuntimeException("Response code was not 200. Detected response was "+responseCode);
			 	}

			 
			 //response = new TeeInputStream( httpConnection.getInputStream(),System.err,false);
			 response =httpResponse.getEntity().getContent();
			 if(this.teeResponse)
				 {
				 stderr().println(queryb);
				 response = new TeeInputStream(response,stderr(),false);
				 }
			 
			
			  
			if(xml_answer)
				{
				return documentBuilder.parse(response);
				}
			else
				{
				inputSource =new StreamSource(response);
				return unmarshaller.unmarshal(inputSource, Opt.class).getValue();
				}
			
			} 
		catch (Exception e)
			{
			throw new IOException(e);
			}
		finally
			{
			CloserUtil.close(response);
			if(httpPost!=null) httpPost.releaseConnection();
			this.lastMillisec = System.currentTimeMillis(); 
			}
		}
	
	private Opt vep(List<VariantContext> contexts) throws IOException
		{
		return (Opt)generic_vep(contexts,false);
		}
	private Document vepxml(List<VariantContext> contexts) throws IOException
		{
		return (Document)generic_vep(contexts,true);
		}
	
	@Override
	public int doWork(final List<String> args) {
	try {
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		this.documentBuilder=dbf.newDocumentBuilder();
		
		final JAXBContext context = JAXBContext.newInstance("org.ensembl.vep");
		this.unmarshaller=context.createUnmarshaller();
		
		TransformerFactory trf=TransformerFactory.newInstance();
		this.xmlSerializer = trf.newTransformer();
		
		/** create http client */
		this.httpClient = HttpClients.createDefault();
		return doVcfToVcf(args, this.outputFile);
		}
	catch(Exception err) {
		LOG.error(err);
		return -1;
	}
	finally 
		{
		this.unmarshaller=null;
		CloserUtil.close(this.httpClient);
		this.httpClient=null;
			
		}
	}
		
	
	@Override
	protected int doVcfToVcf(String inputName, VcfIterator vcfIn, VariantContextWriter out) {
	    try {
		final java.util.Base64.Encoder  base64Encoder=java.util.Base64.getEncoder();
		final SequenceOntologyTree soTree= SequenceOntologyTree.getInstance();
		VCFHeader header=vcfIn.getHeader();
		List<VariantContext> buffer=new ArrayList<>(this.batchSize+1);
		VCFHeader h2= new VCFHeader(header);
		addMetaData(h2);
		
		if(!xmlBase64)
			{
			h2.addMetaDataLine(new VCFInfoHeaderLine(
					TAG,
					VCFHeaderLineCount.UNBOUNDED,
					VCFHeaderLineType.String,
					"VEP Transcript Consequences. Format :(biotype|cdnaStart|cdnaEnd|cdsStart|cdsEnd|geneId|geneSymbol|geneSymbolSource|hgnc|strand|transcript|variantAllele|so_acns)"
					));
			}
		else
			{
			h2.addMetaDataLine(new VCFInfoHeaderLine(
					TAG,
					1,
					VCFHeaderLineType.String,
					"VEP xml answer encoded as base 64"
					));
			}
		
		out.writeHeader(h2);
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
		for(;;)
			{
			VariantContext ctx=null;
			if(vcfIn.hasNext())
				{
				buffer.add((ctx=progress.watch(vcfIn.next())));
				}
			if(ctx==null || buffer.size()>=this.batchSize)
				{
				if(!buffer.isEmpty())
					{
					if(!xmlBase64)
						{
						Opt opt = vep(buffer);
						for(VariantContext ctx2:buffer)
							{
							VariantContextBuilder vcb=new VariantContextBuilder(ctx2);
							final String inputStr = createInputContext(ctx2);
							Data mydata=null;
							for(Data data:opt.getData())
								{
								if(!inputStr.equals(data.getInput())) continue;
								mydata=data;
								break;
								}
							if(mydata==null)
								{
								LOG.info("No Annotation found for "+inputStr);
								out.add(ctx2);
								continue;
								}
							List<String> infoList=new ArrayList<>();
							List<TranscriptConsequences> csql=mydata.getTranscriptConsequences();
							for(int i=0;i< csql.size();++i)
								{
								TranscriptConsequences csq= csql.get(i);
								StringBuilder sb=new StringBuilder();
								sb.append(empty(csq.getBiotype())).append("|").
									append(empty(csq.getCdnaStart())).append("|").
									append(empty(csq.getCdnaEnd())).append("|").
									append(empty(csq.getCdsStart())).append("|").
									append(empty(csq.getCdsEnd())).append("|").
									append(empty(csq.getGeneId())).append("|").
									append(empty(csq.getGeneSymbol())).append("|").
									append(empty(csq.getGeneSymbolSource())).append("|").
									append(empty(csq.getHgncId())).append("|").
									append(empty(csq.getStrand())).append("|").
									append(empty(csq.getTranscriptId())).append("|").
									append(empty(csq.getVariantAllele())).append("|")
										;
								List<String> terms=csq.getConsequenceTerms();
								for(int j=0;j< terms.size();++j)
									{
									if(j>0) sb.append("&");
									SequenceOntologyTree.Term term = soTree.getTermByLabel(terms.get(j));
									if(term==null)
										{
										sb.append(terms.get(j));
										LOG.warn("No SO:Term found for "+terms.get(j));
										}
									else
										{
										sb.append(term.getAcn());
										}
									}
								infoList.add(sb.toString());
								}
							if(!infoList.isEmpty())
								{
								vcb.attribute(TAG, infoList);
								}
							
							out.add(vcb.make());
							}
						}//end of not(XML base 64)
					else
						{
						Document opt = vepxml(buffer);
						Element root= opt.getDocumentElement();
						if(!root.getNodeName().equals("opt"))
							throw new IOException("Bad root node "+root.getNodeName());
						
						for(VariantContext ctx2:buffer)
							{
							String inputStr = createInputContext(ctx2);							
							Document newdom=null;
							
							//loop over <data/>
							for(Node dataNode =root.getFirstChild();
									dataNode!=null;
									dataNode=dataNode.getNextSibling())
								{
								if(dataNode.getNodeType()!=Node.ELEMENT_NODE) continue;
								Attr att = Element.class.cast(dataNode).getAttributeNode("input");
								if(att==null)
									{
									LOG.warn("no @input in <data/>");
									continue;
									}

								if(!att.getValue().equals(inputStr)) continue;
								if(newdom==null)
									{
									newdom = this.documentBuilder.newDocument();
									newdom.appendChild(newdom.createElement("opt"));
									}
								newdom.getDocumentElement().appendChild(newdom.importNode(dataNode, true));
								}
							if(newdom==null)
								{
								LOG.warn("No Annotation found for "+inputStr);
								out.add(ctx2);
								continue;
								}
							StringWriter sw=new StringWriter();
							try {
								this.xmlSerializer.transform(
										new DOMSource(newdom),
										new StreamResult(sw)
										);
								} 
							catch (TransformerException err)
								{
								throw new IOException(err);
								}
							VariantContextBuilder vcb=new VariantContextBuilder(ctx2);
							vcb.attribute(TAG,base64Encoder.encodeToString(sw.toString().getBytes()).
									replaceAll("[\\s=]", ""));
							out.add(vcb.make());
							}
						}//end of XML base 64
					}
				if(ctx==null) break;
				buffer.clear();
				}
			if(out.checkError()) break;
			}
		progress.finish();
		return RETURN_OK;
	    } catch(Exception err)
	    	{
	    	LOG.error(err);
	    	return -1;
	    	}
	
		}
	
	
	public static void main(String[] args) {
		new VcfEnsemblVepRest().instanceMainWithExit(args);
	}
	}
