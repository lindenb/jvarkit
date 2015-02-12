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
* 2015 rewritten using htsjdk
* 2014 creation

*/package com.github.lindenb.jvarkit.tools.vcfgo;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import javax.xml.stream.XMLStreamException;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.go.GoTree;
import com.github.lindenb.jvarkit.util.htsjdk.HtsjdkVersion;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser.SnpEffPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser.VepPrediction;


public class VcfGeneOntology extends AbstractCommandLineProgram
	{
	private String GO="http://archive.geneontology.org/latest-termdb/go_daily-termdb.rdf-xml.gz";
	private String GOA="http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/gene-associations/gene_association.goa_human.gz?rev=HEAD";
	private GoTree goTree=null;
	private Map<String,Set<GoTree.Term>> name2go=null;


	private VcfGeneOntology()
		{
		
		}
	
	@Override
	protected String getOnlineDocUrl()
		{
		return "https://github.com/lindenb/jvarkit/wiki/VCFGeneOntology";
		}
    
    @Override
    public String getProgramDescription()
    	{
    	return "Find the GO terms for VCF annotated with SNPEFF or VEP. ";
    	}

    
	private void readGO() throws IOException
		{
		if(goTree!=null) return;
		this.info("read GO "+GO);
		try
			{
			goTree=GoTree.parse(GO);
			this.info("GO size:"+goTree.size());
			}
		catch(XMLStreamException err)
			{
			throw new IOException(err);
			}
		}
	protected void readGOA() throws IOException
		{
		if(this.name2go!=null) return;
		this.name2go = new HashMap<String, Set<GoTree.Term>>();
		this.info("read GOA "+GOA);
		Pattern tab=Pattern.compile("[\t]");
		BufferedReader in=IOUtils.openURIForBufferedReading(GOA);
		String line;
		while((line=in.readLine())!=null)
			{
			if(line.isEmpty() || line.startsWith("!")) continue;
			String tokens[]=tab.split(line,6);
			if(tokens.length<6) continue;
			
			GoTree.Term term=goTree.getTermByAccession(tokens[4]);
			if(term==null)
				{
				
				continue;
				}
			
			
			Set<GoTree.Term> set=name2go.get(tokens[2]);
			if(set==null)
				{
				set=new HashSet<GoTree.Term>();
				name2go.put(tokens[2],set);
				}
			set.add(term);
			}
		in.close();
		this.info("GOA size:"+name2go.size());
		}
	private String TAG="GOA";
	    
    @Override
    public void printOptions(PrintStream out)
    	{
		out.println("-A (goa input url). Default:"+GOA);
		out.println("-G (go url). Default:"+GO);	
		out.println("-r remove variant if no GO term is found associated to variant");	
		out.println("-T (string) INFO tag. Default:"+TAG);	
		out.println("-C (Go:Term) Add children to the list of go term to be filtered. Can be used multiple times.");	
		out.println("-v inverse filter if -C is used.");	
		out.println("-F (filter) if -C is used, don't remove the variant but set the filter");	
    	super.printOptions(out);
    	}

		
	@Override
	public int doWork(String[] args)
		{
		String filterName=null;
		boolean inverse_filter=false;
		boolean removeIfNoGo=false;
		Set<String> strGoTermToFilter=new HashSet<String>();
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"A:G:rT:C:vF:"))!=-1)
			{
			switch(c)
				{
				case 'T':this.TAG=opt.getOptArg();break;
				case 'A':this.GOA=opt.getOptArg();break;
				case 'G':this.GO=opt.getOptArg();break;
				case 'r':removeIfNoGo=true;break;
				case 'C':strGoTermToFilter.add(opt.getOptArg());break;
				case 'v': inverse_filter=true;break;
				case 'F':filterName=opt.getOptArg();break;

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
		
		
		
		VariantContextWriter w=null;
		VcfIterator in=null;
		try
			{
			readGO();
			readGOA();
			
			Set<GoTree.Term> goTermToFilter=null;
			if(!strGoTermToFilter.isEmpty())
				{
				goTermToFilter=new HashSet<GoTree.Term>(strGoTermToFilter.size());
				for(String acn: strGoTermToFilter)
					{
					GoTree.Term t=this.goTree.getTermByAccession(acn);
					if(t==null)
						{
						error("Cannot find GO acn "+acn);
						return -1;
						}
					goTermToFilter.add(t);
					}
				}
			
			if(opt.getOptInd()==args.length)
				{
				info("reading from stdin.");
				in= VCFUtils.createVcfIteratorStdin();
				}
			else if(opt.getOptInd()+1==args.length)
				{
				String filename=args[opt.getOptInd()];
				info("reading from "+filename);
				in= VCFUtils.createVcfIterator(filename);
				}
			else
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			w=VCFUtils.createVariantContextWriterToStdout();
			VCFHeader header=in.getHeader();
			VCFHeader h2=new VCFHeader(header);
			h2.addMetaDataLine(new VCFInfoHeaderLine(TAG,VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String,"GO terms from GO "+GO+" and GOA="+GOA));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkVersion",HtsjdkVersion.getVersion()));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"HtsJdkHome",HtsjdkVersion.getHome()));
			if(filterName!=null)
				{
				h2.addMetaDataLine( new VCFFilterHeaderLine(
						filterName ,
						"Flag  GO terms "+(inverse_filter?" not descendant of ":"")+" the provided GO terms"
						));
				}
			
			w.writeHeader(h2);
			SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			SnpEffPredictionParser snpEffPredictionParser= new SnpEffPredictionParser(header);
			VepPredictionParser vepPredictionParser = new VepPredictionParser(header);
			
			while(in.hasNext())
				{
				if(System.out.checkError()) break;
				VariantContext ctx=in.next();
				progess.watch(ctx.getChr(),ctx.getStart());
				/* symbols for this variant */
				Set<String> symbols=new HashSet<String>();
				
				/* scan SNPEFF gene */
				for(SnpEffPrediction pred: snpEffPredictionParser.getPredictions(ctx))
					{
					String genName=pred.getGeneName();
					if(genName==null || genName.isEmpty()) continue;
					symbols.add(genName);
					}
				
				/* scan VEP gene */
				for(VepPrediction pred: vepPredictionParser.getPredictions(ctx))
					{
					String genName=pred.getGeneName();
					if(!(genName==null || genName.isEmpty()))
						{
						symbols.add(genName);
						}
					genName=pred.getGene();
					if(!(genName==null || genName.isEmpty()))
						{
						symbols.add(genName);
						}
					genName=pred.getHGNC();
					if(!(genName==null || genName.isEmpty()))
						{
						symbols.add(genName);
						}
					}
				/* only keep known GENES from GOA */
				symbols.retainAll(this.name2go.keySet());
				
				boolean found_child_of_filter=false;
				
				/* ATTS */
				List<String> atts=new ArrayList<String>();
				
				/* loop over symbols */
				for(String symbol:symbols)
					{
					/* go terms associated to this symbol */
					Set<GoTree.Term> t2= this.name2go.get(symbol);
					if(t2==null || t2.isEmpty()) continue;
					
					StringBuilder sb=new StringBuilder(symbol);
					sb.append("|");
					boolean first=true;
					for(GoTree.Term gt:t2)
						{
						/* user gave terms to filter */
						if(!found_child_of_filter && goTermToFilter!=null)
							{
							for(GoTree.Term userTerm :goTermToFilter)
								{
								if(userTerm.hasDescendant(gt.getAcn()))
									{
									found_child_of_filter=true;
									break;
									}
								}
							}
						
						if(!first) sb.append("&");
						sb.append(gt.getAcn());
						first=false;
						}
					atts.add(sb.toString());
					}
				/* no go term was found */
				if(atts.isEmpty())
					{
					if(!removeIfNoGo)
							{
							w.add(ctx);
							}
					continue;
					}
				
				
				VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				
				
				
				
				/* check children of user's terms */
				if(goTermToFilter!=null)
					{
					/* keep if found children*/
					if( ( inverse_filter && found_child_of_filter)||
						(!inverse_filter && !found_child_of_filter))
						{
						/* don't remove, but set filter */
						if(filterName!=null)
							{
							Set<String> filters=new HashSet<String>(ctx.getFilters());
							filters.add(filterName);
							vcb.filters(filters);
							}
						else
							{
							continue;
							}
						}
					}
				
				/* add go terms */
				vcb.attribute(TAG, atts);
				
				w.add(vcb.make());
				}
			progess.finish();
			w.close();
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(w);
			CloserUtil.close(in);
			}
		}
	
		
		public static void main(String[] args)
			{
			new VcfGeneOntology().instanceMainWithExit(args);
			}
		}