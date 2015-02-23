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
import java.io.File;
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
import com.github.lindenb.jvarkit.knime.KnimeApplication;
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

/**
 * 
 * VcfGeneOntology
 *
 */
public class VcfGeneOntology
	extends AbstractCommandLineProgram
	implements KnimeApplication
	{
	private File fileOut=null;
	private int countVariants=0;
	private String GO="http://archive.geneontology.org/latest-termdb/go_daily-termdb.rdf-xml.gz";
	private String GOA="http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/gene-associations/gene_association.goa_human.gz?rev=HEAD";
	private GoTree goTree=null;
	private Map<String,Set<GoTree.Term>> name2go=null;
	private String TAG="GOA";
	private Set<String> strGoTermToFilter=new HashSet<String>();
	private String filterName=null;
	private boolean inverse_filter=false;
	private boolean removeIfNoGo=false;
	private Set<GoTree.Term> goTermToFilter=null;


	/** moved to public for knime */
	public VcfGeneOntology()
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

    public void setGeneOntologyUrl(String gO) {
		this.GO = gO;
		}
    
    public void setGeneOntologyAnnotationUrl(String gOA)
    	{
    	this.GOA = gOA;
		}
    
    public int getVariantCount()
    	{
		return countVariants;
		}
    
    public void setInfoTag(String tAG) {
		this.TAG = tAG;
		}
    
    public void setRemoveIfNoGo(boolean removeIfNoGo) {
		this.removeIfNoGo = removeIfNoGo;
		}
    
    public void setInverseFilter(boolean inverse_filter) {
		this.inverse_filter = inverse_filter;
		}
    
    public void setFilterName(String filterName)
    	{
		this.filterName = filterName;
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
			
			GoTree.Term term= this.goTree.getTermByAccession(tokens[4]);
			if(term==null)
				{
				
				continue;
				}
			
			
			Set<GoTree.Term> set= this.name2go.get(tokens[2]);
			if(set==null)
				{
				set=new HashSet<GoTree.Term>();
				this.name2go.put(tokens[2],set);
				}
			set.add(term);
			}
		in.close();
		this.info("GOA size:"+this.name2go.size());
		}
	
	
	
	
	
	
	@Override
	public int initializeKnime()
		{
		try
			{
			readGO();
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		
		try
			{
			readGOA();
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		
		
		if(!this.strGoTermToFilter.isEmpty())
			{
			goTermToFilter=new HashSet<GoTree.Term>(strGoTermToFilter.size());
			for(String acn: this.strGoTermToFilter)
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
		return 0;
		}

	@Override
	public int executeKnime(List<String> args)
		{
		VcfIterator vcfIn=null;
		try
			{
			if(args.isEmpty())
				{
				vcfIn = VCFUtils.createVcfIteratorStdin();
				}
			else if(args.size()==1)
				{
				vcfIn= VCFUtils.createVcfIterator(args.get(0));
				}
			else
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			this.filterVcfIterator(vcfIn);
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(vcfIn);
			}
		}

	@Override
	public void disposeKnime() {
		
	}

	@Override
	public void checkKnimeCancelled() {
		
	}

	@Override
	public void setOutputFile(File fileOut) {
		this.fileOut=fileOut;
	}
	
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
    
	private void filterVcfIterator(VcfIterator in) throws IOException
		{
		this.countVariants=0;
		VariantContextWriter w = null;
		try {
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
			
			if(this.fileOut==null)
				{
				w = VCFUtils.createVariantContextWriterToStdout();
				}
			else
				{
				info("opening vcf writer to "+this.fileOut);
				w = VCFUtils.createVariantContextWriter(this.fileOut);
				}

			w.writeHeader(h2);
			SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			SnpEffPredictionParser snpEffPredictionParser= new SnpEffPredictionParser(header);
			VepPredictionParser vepPredictionParser = new VepPredictionParser(header);
			
			while(in.hasNext())
				{
				if(System.out.checkError()) break;
				VariantContext ctx=progess.watch(in.next());
				
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
						if(!found_child_of_filter && this.goTermToFilter!=null)
							{
							for(GoTree.Term userTerm : this.goTermToFilter)
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
							this.countVariants++;
							w.add(ctx);
							}
					continue;
					}
				
				
				VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				
				/* check children of user's terms */
				if(this.goTermToFilter!=null)
					{
					/* keep if found children*/
					if( ( this.inverse_filter && found_child_of_filter)||
						(!this.inverse_filter && !found_child_of_filter))
						{
						/* don't remove, but set filter */
						if(this.filterName!=null)
							{
							Set<String> filters=new HashSet<String>(ctx.getFilters());
							filters.add(this.filterName);
							vcb.filters(filters);
							}
						else
							{
							continue;
							}
						}
					}
				
				/* add go terms */
				vcb.attribute(this.TAG, atts);
				w.add(vcb.make());
				this.countVariants++;
				}
			progess.finish();
			w.close();
			w=null;
			}
		finally
			{
			CloserUtil.close(w);
			w=null;
			}
		}

	public void addGoTerm(String term)
		{
		this.strGoTermToFilter.add(term);
		}
		
	@Override
	public int doWork(String[] args)
			{
			com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
			int c;
			while((c=opt.getopt(args,getGetOptDefault()+"A:G:rT:C:vF:"))!=-1)
				{
				switch(c)
					{
					case 'T': this.setInfoTag(opt.getOptArg());break;
					case 'A': this.setGeneOntologyAnnotationUrl(opt.getOptArg());break;
					case 'G': this.setGeneOntologyUrl(opt.getOptArg());break;
					case 'r': this.setRemoveIfNoGo(true);break;
					case 'C': this.addGoTerm(opt.getOptArg());break;
					case 'v': this.setInverseFilter(true);break;
					case 'F': this.setFilterName(opt.getOptArg());break;
					default:
						{
						switch(handleOtherOptions(c, opt, args))
							{
							case EXIT_FAILURE: return -1;
							case EXIT_SUCCESS: return 0;
							default:break;
							}
						}
					}
				}
			
			
			if(this.initializeKnime()!=0)
				{
				return -1;
				}
			List<String> L=new ArrayList<String>();
			for(int i=opt.getOptInd();i<args.length;++i)
				{
				L.add(args[i]);
				}
			return this.executeKnime(L);
			}
	
		
		public static void main(String[] args)
			{
			new VcfGeneOntology().instanceMainWithExit(args);
			}


		}