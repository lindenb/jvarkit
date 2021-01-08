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

*/
package com.github.lindenb.jvarkit.tools.burden;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.stream.IntStream;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.pedigree.Pedigree;
import com.github.lindenb.jvarkit.pedigree.PedigreeParser;
import com.github.lindenb.jvarkit.pedigree.Sample;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;
/**

BEGIN_DOC

## Input 

input is a tab delimited text file with the list of paths to the indexed VCF files.

### List mode

All vcf are stored in a list and we run the all combinations of 'recursion' VCFs.

### Line mode

all the VCFs on the line are combined and tested.

## Example

```
$ find dir -type f -name "*.vcf.gz" > paths.txt
$ java -jar ${JVARKIT_HOME}/dist/combinevcffisher.jar --buffer -p the.ped paths.txt

#genes        n-variants      cases-R cases-A controls-R      controls-A      fisher
G1:G2 13      238     4       90      10      0.0010580863010766208
G1:G3 6       241     1       94      6       0.0029829802538355616
G1:G4 6       241     1       94      6       0.0029829802538355616
G1:G5 6       241     1       94      6       0.0029829802538355616
G2:G3 6       241     1       94      6       0.0029829802538355616
G2:G4 6       241     1       94      6       0.0029829802538355616
G2:G5 6       241     1       94      6       0.0029829802538355616
G3:G4 6       241     1       94      6       0.0029829802538355616
```

END_DOC
*/

@Program(name="combinevcffisher",
description="Combine multiple VCF to perform a 'vertical' fisher test.",
keywords= {"vcf","burden","fisher"},
modificationDate="20200910",
creationDate="20200910"
)
public class CombineVcfFisher extends Launcher {
	private final static int STATUS_AFFECTED = 1;
	private final static int STATUS_UNAFFECTED = 0;
	private final static int STATUS_UNDEFINED = -1;
	private static final Logger LOG = Logger.build(CombineVcfFisher.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-p","--pedigree"},description= PedigreeParser.OPT_DESC,required=true)
	private Path pedigreeFile=null;
	@Parameter(names={"--buffer"},description= "buffer all genes/vcf in memory. Avoid to read each VCF twice or more but consummes mmemory.")
	private boolean do_buffer_vcf = false;
	@Parameter(names={"--recursion"},description= "pool size in 'list mode'. 1, 2 or 3.")
	private int recursion = 2;
	@Parameter(names={"--line"},description= "switch to 'line' mode. Default is 'list' mode. See doc")
	private boolean use_line_mode = false;
	@Parameter(names={"--max-fisher"},description= "max fisher value to be printed.")
	private double max_fisher = 1.0;

	private Pedigree pedigree= null;
	private List<String> genotypeSamples = null;
	private int[] sample_status = null;
	
	
	private class GeneGenotypes {
		private final String geneid;
		private final BitSet genotypes = new BitSet(sample_status.length);
		private int nVariants = 0;
		GeneGenotypes(final String geneid) {
			this.geneid = geneid;
		}
	}
	
	
	private GeneGenotypes loadGeneGenotypes(final Path path) throws IOException {
		try (VCFReader reader = VCFReaderFactory.makeDefault().open(path,false)) {
			final VCFHeader header = reader.getHeader();
			if(!header.hasGenotypingData()) {
				LOG.error("not genotype in "+path);
				return null;
			}
			
			/** this is the first VCF header we ever had */
			if(this.genotypeSamples==null) {
				/** build the list of sample */
				this.genotypeSamples = header.getGenotypeSamples();
				/** build the list of sample status */
				this.sample_status = new int[this.genotypeSamples.size()];
				/** fill list of sample status */
				for(int i=0;i< this.genotypeSamples.size();i++) {
					final String sn = this.genotypeSamples.get(i);
					final Sample sample = this.pedigree.getSampleById(sn);
					int status = STATUS_UNDEFINED;
					if(sample!=null) {
						switch(sample.getStatus()) {
							case affected: status = STATUS_AFFECTED; break;
							case unaffected: status = STATUS_UNAFFECTED; break;
							default : status = STATUS_UNDEFINED; break;
							}
						}
					this.sample_status[i] = status;
					}
				// check we have affected
				if(IntStream.of(this.sample_status).noneMatch(V->V==STATUS_AFFECTED)) {
					LOG.warn("No affected sample in "+path);
					return null;
					}
				// check we have unaffected
				if(IntStream.of(this.sample_status).noneMatch(V->V==STATUS_UNAFFECTED)) {
					LOG.warn("No unaffected sample in "+path);
					return null;
					}
				}
			else if(!this.genotypeSamples.equals( header.getGenotypeSamples())) {
				throw new IOException("samples are not the same with previous VCF and "+path);
				}
			// loop over variants
			final GeneGenotypes geneGenotypes = new GeneGenotypes(IOUtils.getFilenameWithoutCommonSuffixes(path));
			try(CloseableIterator<VariantContext> iter=reader.iterator()) {
				while(iter.hasNext()) {
					final VariantContext ctx = iter.next();
					if(ctx.isFiltered() || !ctx.isVariant()) continue;
					geneGenotypes.nVariants++;
					for(int i=0;i< this.sample_status.length;i++) {
						if(this.sample_status[i]==STATUS_UNDEFINED) continue;
						if(geneGenotypes.genotypes.get(i)) continue;
						final Genotype gt = ctx.getGenotype(i);
						if(gt.isHomRef() || gt.isNoCall()) continue;
						geneGenotypes.genotypes.set(i);
						}
					}
				}
			if(geneGenotypes.nVariants==0 || geneGenotypes.genotypes.nextSetBit(0)==-1) {
				LOG.warn("no valid data for "+path);
				return null;
			}
			return geneGenotypes;
			}
		}
	
	private void runFisher(final PrintWriter pw,final GeneGenotypes ggt) {
		if(ggt.genotypes.nextSetBit(0)==-1) return;
		
		int count_case_sv0 = 0;
		int count_ctrl_sv0 = 0;
		int count_case_sv1 = 0;
		int count_ctrl_sv1 = 0;
	
	
		for(int i=0;i< this.sample_status.length;i++) {
			final int status = this.sample_status[i];
			if(status == STATUS_UNDEFINED) continue;
			final boolean hasVariant = ggt.genotypes.get(i);
			if(!hasVariant) {
				if(status==STATUS_AFFECTED) count_case_sv0++;
				else count_ctrl_sv0++;
			} else // AT_LEAST_ONE_VARIANT 
				{
				if(status==STATUS_AFFECTED) count_case_sv1++;
				else count_ctrl_sv1++;
				}
		}//end of person
	
	
	
	
		final FisherExactTest fisher = FisherExactTest.compute(
				count_case_sv0, count_case_sv1,
				count_ctrl_sv0, count_ctrl_sv1
				);
		
		final double fisherV = fisher.getAsDouble();
		if(fisherV> this.max_fisher) return;
		
		
		pw.print(ggt.geneid);
		pw.print("\t");
		pw.print(ggt.nVariants);
		pw.print("\t");
		pw.print(count_case_sv0);
		pw.print("\t");
		pw.print(count_case_sv1);
		pw.print("\t");
		pw.print(count_ctrl_sv0);
		pw.print("\t");
		pw.print(count_ctrl_sv1);
		pw.print("\t");
		pw.println(fisher.getAsDouble());
		

	}
	
	
	@Override
	public int doWork(final List<String> args) {
		try {
			this.pedigree = new PedigreeParser().parse(this.pedigreeFile);
			
			
			final Map<Path, Optional<GeneGenotypes>> path2genegenotype = new HashMap<>(50_000);
			
			final Function<Path, GeneGenotypes> geneMapper = PATH->{
				
				if(this.do_buffer_vcf && path2genegenotype.containsKey(PATH)) {
					return path2genegenotype.get(PATH).orElse(null);
					}
				final GeneGenotypes value;
				try {
					value = loadGeneGenotypes(PATH);
					}
				catch(final IOException err) {
					throw new RuntimeIOException(err);
					}
				if(this.do_buffer_vcf) {
					path2genegenotype.put(PATH,Optional.ofNullable(value));
					}
				return value;
				};

			 final BinaryOperator<GeneGenotypes>	merger = (A,B) ->{
					final GeneGenotypes ggt = new GeneGenotypes(A.geneid+":"+B.geneid);
					ggt.genotypes.or(A.genotypes);
					ggt.genotypes.or(B.genotypes);
					ggt.nVariants = A.nVariants + B.nVariants;
					return ggt;
				};

				
			
			final List<Path> geneList = new ArrayList<>();
			final String input = oneAndOnlyOneFile(args);
			
			try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				pw.println("#genes\tn-variants\tcases-R\tcases-A\tcontrols-R\tcontrols-A\tfisher");

				
				try(BufferedReader br = super.openBufferedReader(input)) {
					String line;
					while((line=br.readLine())!=null) {
						if(line.startsWith("#") || StringUtils.isBlank(line)) continue;
						final String tokens[] = CharSplitter.TAB.split(line);
						
						if(!use_line_mode) {
							 if(tokens.length!=1) {
									LOG.error("inconsistent number of column in "+line+" exepcted one column but got "+tokens.length);
									return -1;
									}
							geneList.add(Paths.get(tokens[0]));
							}
						else // line mode
							{
							final GeneGenotypes ggt =  Arrays.stream(tokens).
									map(S->Paths.get(S)).
									map(geneMapper).
									filter(X->X!=null).
									reduce(merger).orElse(null);
							if(ggt!=null) runFisher(pw,ggt);
							}
						}
					}
				
				if(!use_line_mode) {
					Collections.sort(geneList,(A,B)->A.getFileName().compareTo(B.getFileName()));
					switch(recursion) {
						case 1: 
							{
							for(int x=0;x< geneList.size();x++) {
								final GeneGenotypes ggtx = geneMapper.apply(geneList.get(x));
								if(ggtx==null) continue;
								runFisher(pw, ggtx);
								}
							break;
							}
						case 2: 
							{
							for(int x=0;x +1< geneList.size();x++) {
								final GeneGenotypes ggtx = geneMapper.apply(geneList.get(x));
								if(ggtx==null) continue;
								for(int y=x+1;y< geneList.size();y++) {
									final GeneGenotypes ggty = geneMapper.apply(geneList.get(y));
									if(ggty==null) continue;
									runFisher(pw, merger.apply(ggtx,ggty));
									}
								}
							break;
							}
						case 3: 
							{
							for(int x=0;x+2< geneList.size();x++) {
								final GeneGenotypes ggtx = geneMapper.apply(geneList.get(x));
								if(ggtx==null) continue;

								for(int y=x+1;y+1< geneList.size();y++) {
									final GeneGenotypes ggty = geneMapper.apply(geneList.get(y));
									if(ggty==null) continue;

									final GeneGenotypes ggtxy = merger.apply(ggtx,ggty);

									
									for(int z=y+1;z< geneList.size();z++) {
										final GeneGenotypes ggtz = geneMapper.apply(geneList.get(z));
										if(ggtz==null) continue;

										
										runFisher(pw, merger.apply(ggtxy,ggtz));
										}
									}
								}
							break;
							}
						default:
							{
							LOG.error("Cannot compute recursion>3 or <1");
							return -1;
							}
						}
					}
				pw.flush();
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	
	public static void main(String[] args) {
		new CombineVcfFisher().instanceMainWithExit(args);
	}

}
