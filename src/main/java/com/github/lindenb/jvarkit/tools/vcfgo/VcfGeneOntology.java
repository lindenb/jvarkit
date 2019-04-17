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
package com.github.lindenb.jvarkit.tools.vcfgo;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import javax.xml.stream.XMLStreamException;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.go.GoTree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import htsjdk.variant.vcf.VCFIterator;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser.SnpEffPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser.VepPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;

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

/**
BEGIN_DOC

## Example


```make
INPUT=${HOME}/test.vcf.gz
VCFGO=java -jar dist/vcfgo.jar -A  gene_association.goa_human.gz -G go_daily-termdb.rdf-xml.gz 
DEPS=dist/vcfgo.jar gene_association.goa_human.gz go_daily-termdb.rdf-xml.gz 

all: $(addsuffix .vcf,$(addprefix test,1 2 3 4 5 6 7))
	$(foreach F,$^,echo -e "Output of $F:"  && cat $F; )

test1.vcf: ${DEPS}
	${VCFGO} ${INPUT} |  grep -v "#" | grep GO | cut -f 7,8 | head -n 1 > $@

test2.vcf: ${DEPS}
	${VCFGO}  -T MYGOATAG ${INPUT} |  grep -v "#" | grep MYGOATAG | cut -f 7,8 | head -n 1 > $@
	
test3.vcf: ${DEPS}
	${VCFGO}  -C GO:0007283 ${INPUT} |  grep -v "#" | grep GO  | cut -f 7,8 | head -n 1 > $@
	
test4.vcf: ${DEPS}
	${VCFGO} -C GO:0007283 -v  ${INPUT}|  grep -v "#" |  grep GO | cut -f 7,8 | head -n 1 > $@

test5.vcf: ${DEPS}
	${VCFGO}  -C GO:0007283 -F GOFILTER ${INPUT}| grep -v "#" |   grep GOFILTER | cut -f 7,8 | head -n1 > $@
	
test6.vcf: ${DEPS}
	${VCFGO} -C GO:0007283 -F GOFILTER -v  ${INPUT}| grep -v "#" |   grep GOFILTER | cut -f 7,8 | head -n1 > $@	
	
test7.vcf: ${DEPS}
	${VCFGO} -r ${INPUT} | grep -v "#" | cut -f 7,8 | head -n1 > $@
		
	
gene_association.goa_human.gz :
	curl -o $@ "http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/gene-associations/gene_association.goa_human.gz?rev=HEAD"

go_daily-termdb.rdf-xml.gz :
	curl -o $@ "http://archive.geneontology.org/latest-termdb/$@"

dist/vcfgo.jar:
	$(MAKE) vcfgo

```
output of GNU make:

Output of test1.vcf:

```
.	(...)GOA=AURKAIP1|GO:0070124&GO:0005739&GO:0005515&GO:0070125&GO:0045839&GO:0070126&GO:0032543&GO:0006996&GO:0005743&GO:0005654&GO:0005634&GO:0045862&GO:0043231;MQ=39	GT:PL:DP:GQ	1/1:35,3,0:1:4
```

Output of test2.vcf:

```
.	(...)MYGOATAG=AURKAIP1|GO:0070124&GO:0005739&GO:0005515&GO:0070125&GO:0045839&GO:0070126&GO:0032543&GO:0006996&GO:0005743&GO:0005654&GO:0005634&GO:0045862&GO:0043231
```

Output of test3.vcf:

```
.	(...)GOA=GJA9|GO:0007154&GO:0005922&GO:0016021,MYCBP|GO:0005739&GO:0005515&GO:0006355&GO:0005813&GO:0005737&GO:0003713&GO:0005634&GO:0007283&GO:0006351;MQ=46
```

Output of test4.vcf:
````
.	(...)GOA=AURKAIP1|GO:0070124&GO:0005739&GO:0005515&GO:0070125&GO:0045839&GO:0070126&GO:0032543&GO:0006996&GO:0005743&GO:0005654&GO:0005634&GO:0045862&GO:0043231
```

Output of test5.vcf:

```
GOFILTER	(...)GOA=AURKAIP1|GO:0070124&GO:0005739&GO:0005515&GO:0070125&GO:0045839&GO:0070126&GO:0032543&GO:0006996&GO:0005743&GO:0005654&GO:0005634&GO:0045862&GO:0043231
```

Output of test6.vcf:

```
GOFILTER	(...)GOA=GJA9|GO:0007154&GO:0005922&GO:0016021,MYCBP|GO:0005739&GO:0005515&GO:0006355&GO:0005813&GO:0005737&GO:0003713&GO:0005634&GO:0007283&GO:0006351
```

Output of test7.vcf:

```
.	(...)GOA=AURKAIP1|GO:0070124&GO:0005739&GO:0005515&GO:0070125&GO:0045839&GO:0070126&GO:0032543&GO:0006996&GO:0005743&GO:0005654&GO:0005634&GO:0045862&GO:0043231
```
END_DOC
 *
 */
@Program(
		name="vcfgo",
		description="Find the GO terms for VCF annotated with SNPEFF or VEP",
		deprecatedMsg="do not use this",
		keywords={"vcf","go"}
		)
public class VcfGeneOntology
	extends Launcher
	{
	 private static Logger LOG=Logger.build(VcfGeneOntology.class).make(); 
		
	 @Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	 private File outputFile = null;

	@Parameter(names="-G",description="(go  url)",required=true)
	private String GO="http://archive.geneontology.org/latest-termdb/go_daily-termdb.rdf-xml.gz";
	@Parameter(names="-A",description="(goa input url)",required=true)
	private String GOA="http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/gene-associations/gene_association.goa_human.gz?rev=HEAD";
	private GoTree goTree=null;
	private Map<String,Set<GoTree.Term>> name2go=null;
	@Parameter(names="-T",description="INFO tag.")
	private String TAG="GOA";
	@Parameter(names="-C",description="(Go:Term) Add children to the list of go term to be filtered. Can be used multiple times.")
	private Set<String> strGoTermToFilter=new HashSet<String>();
	@Parameter(names="-F",description=" if -C is used, don't remove the variant but set the filter")
	private String filterName=null;
	@Parameter(names="-v",description="inverse filter if -C is used")
	private boolean inverse_filter=false;
	@Parameter(names="-r",description="remove variant if no GO term is found associated to variant")
	private boolean removeIfNoGo=false;

	private Set<GoTree.Term> goTermToFilter=null;

	

	/** moved to public for knime */
	public VcfGeneOntology()
		{
		
		}
	

    
	private void readGO() throws IOException
		{
		if(goTree!=null) return;
		LOG.info("read GO "+GO);
		try
			{
			goTree=GoTree.parse(GO);
			LOG.info("GO size:"+goTree.size());
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
		LOG.info("read GOA "+GOA);
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
		LOG.info("GOA size:"+this.name2go.size());
		}
	
	
	
	
	
	
	private int initializeThings()
		{
		try
			{
			readGO();
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		
		try
			{
			readGOA();
			}
		catch(Exception err)
			{
			LOG.error(err);
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
					LOG.error("Cannot find GO acn "+acn);
					return -1;
					}
				goTermToFilter.add(t);
				}
			}
		return 0;
		}

	public int executeThings(List<String> args)
		{
		VCFIterator vcfIn=null;
		try
			{
			vcfIn = super.openVCFIterator(oneFileOrNull(args));
			this.filterVCFIterator(vcfIn);
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(vcfIn);
			}
		}

    
	private void filterVCFIterator(final VCFIterator in) throws IOException
		{
		VariantContextWriter w = null;
		try {
			VCFHeader header=in.getHeader();
			VCFHeader h2=new VCFHeader(header);
			h2.addMetaDataLine(new VCFInfoHeaderLine(TAG,VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String,"GO terms from GO "+GO+" and GOA="+GOA));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			JVarkitVersion.getInstance().addMetaData(getClass().getSimpleName(), h2);
			
			if(filterName!=null)
				{
				h2.addMetaDataLine( new VCFFilterHeaderLine(
						filterName ,
						"Flag  GO terms "+(inverse_filter?" not descendant of ":"")+" the provided GO terms"
						));
				}
			
			w = super.openVariantContextWriter(outputFile);

			w.writeHeader(h2);
			final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			final SnpEffPredictionParser snpEffPredictionParser= new SnpEffPredictionParserFactory().header(header).get();
			final VepPredictionParser vepPredictionParser = new VepPredictionParserFactory().header(header).get();
			
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
				for(final String symbol:symbols)
					{
					/* go terms associated to this symbol */
					final Set<GoTree.Term> t2= this.name2go.get(symbol);
					if(t2==null || t2.isEmpty()) continue;
					
					final StringBuilder sb=new StringBuilder(symbol);
					sb.append("|");
					boolean first=true;
					for(final GoTree.Term gt:t2)
						{
						/* user gave terms to filter */
						if(!found_child_of_filter && this.goTermToFilter!=null)
							{
							for(final GoTree.Term userTerm : this.goTermToFilter)
								{
								if(gt.isDescendantOf(userTerm))
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
		public int doWork(final List<String> args) {			
			if(this.initializeThings()!=0)
				{
				return -1;
				}
			return this.executeThings(args);
			}
	
		
		public static void main(String[] args)
			{
			new VcfGeneOntology().instanceMainWithExit(args);
			}


		}
