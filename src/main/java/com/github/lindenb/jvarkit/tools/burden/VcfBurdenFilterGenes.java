/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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

*/
package com.github.lindenb.jvarkit.tools.burden;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.BcfToolsPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.BcfToolsPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;

import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
/**

BEGIN_DOC

### Example

```
echo "IL2" > genes.txt
echo "NOCTH2" >>  genes.txt
gunzip -c input.vcf.gz |\
	java -jar dit/vcfburdenfiltergenes.jar -g genes.txt
```

### Example

```
$ wget -O - -q "https://github.com/immune-health/antigen.garnish/raw/f0453336a4859c83d27640c286d3960c1672f164/inst/extdata/testdata/antigen.garnish_hg19anno_example.vcf"  |\
	grep -w 216371854  | cut -f 8 | tr ";" "\n"   | grep ^ANN= | cut -c5- | tr "," "\n"
T|missense_variant|MODERATE|USH2A|USH2A|transcript|NM_206933.2|protein_coding|18/72|c.3884G>A|p.Arg1295Gln|4271/18883|3884/15609|1295/5202||
T|missense_variant|MODERATE|USH2A|USH2A|transcript|NM_007123.5|protein_coding|18/21|c.3884G>A|p.Arg1295Gln|4271/6316|3884/4641|1295/1546||
```

```
$ wget -O - -q "https://github.com/immune-health/antigen.garnish/raw/f0453336a4859c83d27640c286d3960c1672f164/inst/extdata/testdata/antigen.garnish_hg19anno_example.vcf"  |\
	sed 's/PASS\t[A-Z0-9]*;/PASS\t/' |\ ## the vcf above is malformed, quick hack
	java -jar dist/vcfburdenfiltergenes.jar -a "NM_206933.2" |\
	grep -w 216371854  | cut -f 8 | tr ";" "\n"   | grep ^ANN= | cut -c5- | tr "," "\n"
T|missense_variant|MODERATE|USH2A|USH2A|transcript|NM_206933.2|protein_coding|18/72|c.3884G>A|p.Arg1295Gln|4271/18883|3884/15609|1295/5202||
```


END_DOC
*/
@Program(
		name="vcfburdenfiltergenes",
		description="Filter VEP/SnpEff Output from a list of genes.",
		keywords={"gene","vcf","vep","snpeff"},
		biostars=353011,
		modificationDate="20200924",
		creationDate="20160322"
		)
public class VcfBurdenFilterGenes
	extends OnePassVcfLauncher
	{
	private static final Logger LOG = Logger.build(VcfBurdenFilterGenes.class).make();


	@Parameter(names={"-g","--genes"},description="Gene/transcript file: one name per line")
	private Path geneFile = null;
	@Parameter(names={"-a","--add"},description="[20180627] Gene Names: Add this gene, multiple separated by comma,spaces,semicolon")
	private String geneStr="";
	@Parameter(names={"-filter","--filter"},description="If empty: remove the variants from the VCF. If not empty, add a token in the column FILTER.")
	private String filterTag = "";

	
	private final Set<String> geneNames= new HashSet<>();
	
	
	
	public VcfBurdenFilterGenes()
		{
		}
	
	
	
	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VCFIterator in,
			final VariantContextWriter out
			) 
		{
		final VCFHeader header=in.getHeader();		
		try {
			final VCFHeader h2=addMetaData(new VCFHeader(header));

			final VCFFilterHeaderLine filterControlsHeader;
			if(!StringUtil.isBlank(this.filterTag))
				{
				filterControlsHeader = new VCFFilterHeaderLine(
						this.filterTag.trim(),
						"Genes not in list "+this.geneFile
						);
				h2.addMetaDataLine(filterControlsHeader);
				}
			else
				{
				filterControlsHeader = null;
				}
			final List<String> lookColumns=Arrays.asList(
					"CCDS",
					"Feature",
					"ENSP",
					"Gene",
					"HGNC",
					"HGNC_ID",					
					"SYMBOL",					
					"RefSeq"					
					);
			final VepPredictionParser vepParser = new VepPredictionParserFactory(header).get();
			final AnnPredictionParser annParser = new AnnPredictionParserFactory(header).get();
			final BcfToolsPredictionParser bcftoolsParser = new BcfToolsPredictionParserFactory(header).get();
			JVarkitVersion.getInstance().addMetaData(this, h2);
			out.writeHeader(h2);
			while(in.hasNext())
				{
				final VariantContext ctx = in.next();
				boolean keep=false;
				final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
				
				
				//not just set FILTER ?
				if(filterControlsHeader==null) {
					vcb.rmAttribute(vepParser.getTag());
					vcb.rmAttribute(annParser.getTag());
				}
				
				final List<String> newVepList=new ArrayList<>();
				for(final String predStr: ctx.getAttributeAsStringList(vepParser.getTag(),""))
					{
					final VepPredictionParser.VepPrediction pred = vepParser.parseOnePrediction(ctx,predStr);
					for(final String col:lookColumns) {
						final String token = pred.getByCol(col);
						if(!StringUtil.isBlank(token) && this.geneNames.contains(token))
							{
							newVepList.add(predStr);
							keep=true;
							break;// break lookColumns
							}
						}
					}
				
				final List<String> newEffList=new ArrayList<>();
				for(final String predStr: ctx.getAttributeAsStringList(annParser.getTag(),"")) {
					final AnnPredictionParser.AnnPrediction pred = annParser.parseOnePrediction(predStr);
					String token = pred.getGeneName();
					if(!StringUtil.isBlank(token) && this.geneNames.contains(token))
						{
						newEffList.add(predStr);
						keep=true;
						continue;
						}
					token = pred.getGeneId();
					if(!StringUtil.isBlank(token) && this.geneNames.contains(token))
						{
						newEffList.add(predStr);
						keep=true;
						continue;
						}
					token = pred.getFeatureId();
					if(!StringUtil.isBlank(token) && this.geneNames.contains(token))
						{
						newEffList.add(predStr);
						keep=true;
						continue;
						}
					}
				
				final List<String> newBcfList=new ArrayList<>();
				for(final String predStr: ctx.getAttributeAsStringList(bcftoolsParser.getTag(), "")) {
					final BcfToolsPredictionParser.BcfToolsPrediction pred = bcftoolsParser.parseOnePrediction(ctx,predStr);
					String token = pred.getGeneName();
					if(!StringUtil.isBlank(token) && this.geneNames.contains(token))
						{
						newBcfList.add(predStr);
						keep=true;
						continue;
						}
					token = pred.getTranscript();
					if(!StringUtil.isBlank(token) && this.geneNames.contains(token))
						{
						newBcfList.add(predStr);
						keep=true;
						continue;
						}
					}
				
				
				//not just set FILTER ?
				if(filterControlsHeader==null) {
					if(!newVepList.isEmpty()) vcb.attribute(vepParser.getTag(),newVepList);
					if(!newEffList.isEmpty()) vcb.attribute(annParser.getTag(),newEffList);
					if(!newBcfList.isEmpty()) vcb.attribute(bcftoolsParser.getTag(),newBcfList);
					}
				
				
				if(filterControlsHeader!=null)
					{
					if(!keep)
						{
						vcb.filter(filterControlsHeader.getID());
						}
					else if(!ctx.isFiltered())
						{
						vcb.passFilters();
						}
					out.add(vcb.make());
					}
				else
					{
					if(keep) out.add(vcb.make());
					}
				}
			return 0;
			} catch(final Throwable err) {
				LOG.error(err);
				return -1;
			}
		}
	
	 
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected int beforeVcf() {
		if(StringUtil.isBlank(this.geneStr) && (this.geneFile==null || !Files.exists(this.geneFile))) {
			LOG.error("Undefined gene file option.");
			return -1;
			}
		
		
		this.geneNames.clear();
		if(this.geneFile!=null)
			{
			try {
				this.geneNames.addAll(Files.readAllLines(this.geneFile));
				}
			catch(IOException err) {
				LOG.error(err);
				return -1;
				}
			}
		Arrays.stream(this.geneStr.split("[ ;\t\n,]+")).
			filter(gs->!StringUtil.isBlank(gs)).
			forEach(gs->this.geneNames.add(gs));
			
		
		this.geneNames.remove(".");
		this.geneNames.remove("");
		LOG.info("number of genes : "+geneNames.size());

		return 0;
		}
	
	
	
	public static void main(final String[] args)
		{
		new VcfBurdenFilterGenes().instanceMainWithExit(args);
		}
	}
