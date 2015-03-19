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

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import javax.xml.transform.stream.StreamSource;

import org.ensembl.vep.*;

import com.github.lindenb.jvarkit.io.TeeInputStream;
import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter3;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class VcfEnsemblVepRest extends AbstractVCFFilter3
	{
	@SuppressWarnings("unused")
	private static final ObjectFactory _fool_javac=null;
	
	private String server = "http://grch37.rest.ensembl.org";
	private String extension = "/vep/homo_sapiens/region";

	private Unmarshaller unmarshaller=null;
	
	@Override
	public String getProgramDescription() {
		return "";
		}
	
	public void setServer(String server) {
		this.server = server;
		}
	
	public void setExtension(String extension) {
		this.extension = extension;
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -o (fileout) output. Default: stdout");
		out.println(" -s (server) REST server . Default: "+this.server);
		out.println(" -e (path extension)  Default: stdout"+this.extension);
		super.printOptions(out);
		}
	
	private static String createInputContext(VariantContext ctx)
		{
		StringBuilder sb=new StringBuilder();
		sb.append(ctx.getChr()).
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
	
	private Opt vep(List<VariantContext> contexts) throws IOException
		{
		info("Running VEP "+contexts.size());
		OutputStream wr=null;
		URLConnection urlConnection=null;
		HttpURLConnection  httpConnection=null;
		InputStream response =null;
		javax.xml.transform.Source inputSource=null;
		try {
			 URL url = new URL(this.server + this.extension);
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
			 System.err.println(queryb);
			 urlConnection = url.openConnection();
			 httpConnection = (HttpURLConnection)urlConnection;
			 httpConnection.setRequestMethod("POST");
			 httpConnection.setRequestProperty("Content-Type", "application/json");
			 httpConnection.setRequestProperty("Accept", "text/xml");
			 httpConnection.setRequestProperty("Content-Length", Integer.toString(postBody.length));
			 httpConnection.setUseCaches(false);
			 httpConnection.setDoInput(true);
			 httpConnection.setDoOutput(true);
			  
			 wr = httpConnection.getOutputStream();
			 wr.write(postBody);
			 wr.flush();
			 wr.close();
			 wr=null;
			  
			 response = new TeeInputStream( httpConnection.getInputStream(),System.err,false);
			 int responseCode = httpConnection.getResponseCode();
			  
			 if(responseCode != 200)
			 	{
				throw new RuntimeException("Response code was not 200. Detected response was "+responseCode);
			 	}
			
			  
			
			inputSource =new StreamSource(response);
			Opt opt=	unmarshaller.unmarshal(inputSource, Opt.class).getValue();
			
			response.close();
			response=null;
			httpConnection.disconnect();
			httpConnection=null;
			
			return opt;
			} 
		catch (JAXBException e)
			{
			throw new IOException(e);
			}
		finally
			{
			CloserUtil.close(wr);
			CloserUtil.close(response);
			if(httpConnection!=null) httpConnection.disconnect();
			}
		}
	
	@Override
	protected void doWork(String inputSource, VcfIterator vcfIn,
			VariantContextWriter out) throws IOException
		{
		final int BUFFER_MAX_SIZE=10;
		final SequenceOntologyTree soTree= SequenceOntologyTree.getInstance();
		VCFHeader header=vcfIn.getHeader();
		List<VariantContext> buffer=new ArrayList<>(BUFFER_MAX_SIZE);
		VCFHeader h2= new VCFHeader(header);
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
		h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));
		h2.addMetaDataLine(new VCFInfoHeaderLine(
				"VEPTRCSQ",
				VCFHeaderLineCount.UNBOUNDED,
				VCFHeaderLineType.String,
				"VEP Transcript Consequences. Format :(biotype|cdnaStart|cdnaEnd|cdsStart|cdsEnd|geneId|geneSymbol|geneSymbolSource|hgnc|strand|transcript|variantAllele|so_acns)"
				));

		
		out.writeHeader(h2);
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
		for(;;)
			{
			VariantContext ctx=null;
			if(vcfIn.hasNext())
				{
				buffer.add((ctx=progress.watch(vcfIn.next())));
				}
			if(ctx==null || buffer.size()>=BUFFER_MAX_SIZE)
				{
				if(!buffer.isEmpty())
					{
					Opt opt = vep(buffer);
					for(VariantContext ctx2:buffer)
						{
						VariantContextBuilder vcb=new VariantContextBuilder(ctx2);
						String inputStr = createInputContext(ctx2);
						Data mydata=null;
						for(Data data:opt.getData())
							{
							if(!inputStr.equals(data.getInput())) continue;
							mydata=data;
							break;
							}
						if(mydata==null)
							{
							info("No Annotation found for "+inputStr);
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
									warning("No SO:Term found for "+terms.get(j));
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
							vcb.attribute("VEPTRCSQ", infoList);
							}
						
						out.add(vcb.make());
						}
					}
				if(ctx==null) break;
				buffer.clear();
				}
			}
		progress.finish();
		}
	
	@Override
	public int initializeKnime() {
		JAXBContext context;
		try {
			context = JAXBContext.newInstance("org.ensembl.vep");
			this.unmarshaller=context.createUnmarshaller();
		} catch (JAXBException e) {
			error(e);
			return -1;
		}
		return super.initializeKnime();
		}
	
	@Override
	public void disposeKnime() {
		this.unmarshaller=null;
		super.disposeKnime();
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"o:s:e:"))!=-1)
			{
			switch(c)
				{
				case 's': setServer(opt.getOptArg()); break;
				case 'e': setExtension(opt.getOptArg()); break;
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
		return this.mainWork(opt.getOptInd(), args);
		}
	public static void main(String[] args) {
		new VcfEnsemblVepRest().instanceMainWithExit(args);
	}
	}
