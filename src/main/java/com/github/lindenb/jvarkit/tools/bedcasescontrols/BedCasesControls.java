/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.bedcasescontrols;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.stats.FisherCasesControls;
import com.github.lindenb.jvarkit.pedigree.CasesControls;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.JexlVariantPredicate;
import com.github.lindenb.jvarkit.variant.vcf.BufferedVCFReader;
import com.github.lindenb.jvarkit.variant.vcf.PerContigVcfReader;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;

/**
BEGIN_DOC

## Example


END_DOC

 */
@Program(
		name="bedcasescontrols",
		description="Count number of cases/controls in from a BED and a VCF file",
		keywords={"bed","burden","vcf","case","control"},
		creationDate="20240603",
		modificationDate="20240603",
		jvarkit_amalgamion =  true
		)
public class BedCasesControls
	extends Launcher
	{
	private static final Logger LOG = Logger.of(BedCasesControls.class);
	
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile= null;
	@Parameter(names={"--vcf","--variants"},description="Indexed VCF file")
	private Path vcfPath= null;
	@Parameter(names={"--buffer"},description=BufferedVCFReader.OPT_BUFFER_DESC)
	private int bufferSize=10_000;
	@Parameter(names={"--include"},description=JexlVariantPredicate.PARAMETER_DESCRIPTION,splitter = NoSplitter.class,converter = JexlVariantPredicate.Converter.class)
	private Predicate<VariantContext> jexlVariantPredicate= JexlVariantPredicate.create("");
	@Parameter(names={"--bed12"},description="input is BED12")
	private boolean is_bed12_flag=false;
	@ParametersDelegate
	private CasesControls casesControls = new CasesControls();
	
	
	
	

	private boolean accept(final VariantContext ctx) {
		return jexlVariantPredicate.test(ctx);
	}
	

	private void scan(
			final BufferedReader br,
			final PrintWriter w,
			final VCFReader vcfReader,
			final CasesControls casesControls
			) throws IOException {
		final BEDCodec bed12codec=new BEDCodec();
		final BedLineCodec bedcodec=new BedLineCodec();
		String line;
		final Set<String> all = casesControls.getAll();
		while((line=br.readLine())!=null) {
			if(StringUtils.isBlank(line) || BedLine.isBedHeader(line)) {
				continue;
				}
			final List<Locatable> intervals;
			if(this.is_bed12_flag) {
				final BEDFeature rec = bed12codec.decode(line);
				if(rec==null) continue;
				intervals  = rec.getExons().stream().map(EX->new SimpleInterval(rec.getContig(),EX.getCdStart(),EX.getCdEnd())).collect(Collectors.toList());					
				}
			else
				{
				final BedLine rec=bedcodec.decode(line);
				if(rec==null) continue;
				intervals = Collections.singletonList(rec);
				}
			int n_variants=0;
			final FisherCasesControls fisherCasesControls=new FisherCasesControls(this.casesControls);
			for(Locatable loc:intervals) {
				try(CloseableIterator<VariantContext> iter= vcfReader.query(loc)) {
					while(iter.hasNext()) {
						final VariantContext ctx= iter.next();
						if(!accept(ctx)) continue;
						n_variants++;
						for(final String sn:all) {
							final Genotype g = ctx.getGenotype(sn);
							fisherCasesControls.accept(g);
							}
						}
						
					}
				}
			w.print(line);
			w.print("\t");
			w.print(fisherCasesControls.getCasesAltCount());
			w.print("\t");
			w.print(fisherCasesControls.getCasesRefCount());
			w.print("\t");
			w.print(fisherCasesControls.getControlsAltCount());
			w.print("\t");
			w.print(fisherCasesControls.getControlsRefCount());
			w.print("\t");
			w.print(n_variants);
			w.print("\t");
			w.print(fisherCasesControls.getAsDouble());
			w.println();
			}
		}
	
			
	private VCFReader openVcf(Path p) throws IOException {
		if(p.getFileName().toString().endsWith(".list")) {
			return new PerContigVcfReader(p);
		}
		else
		{
			return new VCFFileReader(this.vcfPath, true);	
		}
	}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			List<Path> bedPaths = IOUtils.unrollPaths(args);
			try(VCFReader vcfFileReader = openVcf(this.vcfPath)) {
				final VCFHeader header = vcfFileReader.getHeader();
				this.casesControls.load().retain(header);
				this.casesControls.checkHaveCasesControls();
				try(VCFReader bufferedVCF = new BufferedVCFReader(vcfFileReader, this.bufferSize)) {
					try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
					if(bedPaths.isEmpty()) {
							try(BufferedReader br=IOUtils.openStdinForBufferedReader()) {
								scan(br,out,bufferedVCF,this.casesControls);
								}
							}
						else
							{
							for(Path p:bedPaths ) {
								try(BufferedReader br=IOUtils.openPathForBufferedReading(p)) {
									scan(br,out,bufferedVCF,this.casesControls);
									}
								}
							}
						out.flush();
						}
					}
				}
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		}
	

	public static void main(final String[] args)
		{
		new BedCasesControls().instanceMainWithExit(args);
		}
	}
