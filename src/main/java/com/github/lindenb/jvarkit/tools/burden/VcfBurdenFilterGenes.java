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


History:

*/
package com.github.lindenb.jvarkit.tools.burden;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

import htsjdk.variant.vcf.VCFIterator;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
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



## History

  * 20181205 : snpeff scan transcriptID
  * 20180617 : for SNpEFF, now looks into GeneName OR GeneId (was only GeneName)

END_DOC
*/
@Program(
		name="vcfburdenfiltergenes",
		description="Filter VEP/SnpEff Output from a list of genes.",
		keywords={"gene","vcf","vep","snpeff"},
		biostars=353011
		)
public class VcfBurdenFilterGenes
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfBurdenFilterGenes.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-g","--genes"},description="Gene/transcript file: one name per line")
	private File geneFile = null;
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
			final ProgressFactory.Watcher<VariantContext> progess=ProgressFactory.newInstance().dictionary(header).logger(LOG).build();
			JVarkitVersion.getInstance().addMetaData(this, h2);
			out.writeHeader(h2);
			while(in.hasNext() &&  !out.checkError())
				{
				final VariantContext ctx = progess.apply(in.next());
				boolean keep=false;
				final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
				
				
				//not just set FILTER ?
				if(filterControlsHeader==null) {
					vcb.rmAttribute(vepParser.getTag());
					vcb.rmAttribute(annParser.getTag());
				}
				
				final List<String> newVepList=new ArrayList<>();
				for(final String predStr: ctx.getAttributeAsList(vepParser.getTag()).stream().map(O->String.class.cast(O)).collect(Collectors.toList()))
					{
					final VepPredictionParser.VepPrediction pred = vepParser.parseOnePrediction(ctx,predStr);
					for(final String col:lookColumns) {
						final String token = pred.getByCol(col);
						if(!StringUtil.isBlank(token) && this.geneNames.contains(token))
							{
							newVepList.add(predStr);
							keep=true;
							break;
							}
						}
					}
				
				final List<String> newEffList=new ArrayList<>();
				for(final String predStr: ctx.getAttributeAsList(annParser.getTag()).stream().map(O->String.class.cast(O)).collect(Collectors.toList())) {
					final AnnPredictionParser.AnnPrediction pred = annParser.parseOnePrediction(predStr);
					String token = pred.getGeneName();
					if(!StringUtil.isBlank(token) && this.geneNames.contains(token))
						{
						newEffList.add(predStr);
						keep=true;
						break;
						}
					token = pred.getGeneId();
					if(!StringUtil.isBlank(token) && this.geneNames.contains(token))
						{
						newEffList.add(predStr);
						keep=true;
						break;
						}
					token = pred.getFeatureId();
					if(!StringUtil.isBlank(token) && this.geneNames.contains(token))
						{
						newEffList.add(predStr);
						keep=true;
						break;
						}
					}
				
				
				
				//not just set FILTER ?
				if(filterControlsHeader==null) {
					if(!newVepList.isEmpty()) vcb.attribute(vepParser.getTag(),newVepList);
					if(!newEffList.isEmpty()) vcb.attribute(annParser.getTag(),newEffList);
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
			progess.close();
			return RETURN_OK;
			} catch(final Exception err) {
				LOG.error(err);
				return -1;
			} finally {
				CloserUtil.close(in);
			}
		}
	
	 	
	@Override
	public int doWork(final List<String> args) {
		if(StringUtil.isBlank(this.geneStr) && (this.geneFile==null || !this.geneFile.exists())) {
			LOG.error("Undefined gene file option.");
			return -1;
			}
			
		try {
			this.geneNames.clear();
			if(this.geneFile!=null)
				{
				this.geneNames.addAll(Files.readAllLines(this.geneFile.toPath()));
				}
			
			for(final String gs: this.geneStr.split("[ ;\t\n,]+"))
				{
				if(StringUtil.isBlank(gs)) continue;
				this.geneNames.add(gs);
				}
			
			geneNames.remove(".");
			geneNames.remove("");
			LOG.info("number of genes : "+geneNames.size());
			return doVcfToVcf(args,this.outputFile);
		} catch (final Exception err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(final String[] args)
		{
		new VcfBurdenFilterGenes().instanceMainWithExit(args);
		}
	}
