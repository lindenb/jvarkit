package com.github.lindenb.jvarkit.tools.regenie;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.OptionalDouble;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

/** 
BEGIN_DOC



END_DOC
*/
@Program(name="regeniebedannot",
description="Create annotation files for regenie using sliding annotations",
keywords={"vcf","regenie","burden"},
creationDate="20250311",
modificationDate="20250311",
generate_doc = false
)
public class RegenieBedAnnot extends AbstractRegenieAnnot {
	private static final Logger LOG = Logger.build(RegenieBedAnnot.class).make();
	@Parameter(names = {"-B","--bed"}, description = "custom bed file chrom/start/end/name[/score]",required = true)
	private Path userBed;
	@Parameter(names = {"-A","--annotation"}, description = "value for annotation field")
	private String annotation_value="";

	
	private final IntervalTreeMap<UserBed> interval2userbed = new IntervalTreeMap<>();

	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	private static class UserBed {
		String name;
		double score;
	}
	
	@Override
	protected VCFHeader initVcfHeader(VCFHeader h) {
		if(StringUtils.isBlank(this.annotation_value)) {
			this.annotation_value=IOUtils.getFilenameWithoutCommonSuffixes(this.userBed).replaceAll("[^a-zA-Z0-9]+","_");
			}
		try(BufferedReader br=IOUtils.openPathForBufferedReading(this.userBed)) {
			final BedLineCodec bc = new  BedLineCodec();
			String line;
			while((line=br.readLine())!=null) {
				if(BedLine.isBedHeader(line)) continue;
				final BedLine rec = bc.decode(line);
				if(rec==null) continue;
				final String ctg = fixContig(rec.getContig());
				if(StringUtils.isBlank(rec.get(3))) throw new IOException("empty title in bed line "+line);
				final Interval r= new Interval(ctg, rec.getStart(), rec.getEnd(), false, rec.get(3));
				final UserBed ub=new UserBed();
				ub.name = r.getName();
				ub.score = Double.parseDouble(rec.getOrDefault(4,"1.0"));
				this.interval2userbed.put(r, ub);
				}
			}
		catch(IOException err) {
			throw new RuntimeIOException(err);
			}
		return super.initVcfHeader(h);
		}
		
	@Override
	protected void dump(PrintWriter w,final VariantContext ctx) throws Exception {	
		for(UserBed ub: this.interval2userbed.getOverlapping(ctx) ) {
			final Variation v = new Variation();
			v.contig = fixContig(ctx.getContig());
			v.pos = ctx.getStart();
			v.id = makeID(ctx);
			v.gene = ub.name;
			v.prediction = this.annotation_value;
			v.score = OptionalDouble.of(ub.score);
			v.cadd = getCaddScore(ctx);
			print(w,v);
			}
		}

	
	

	public static void main(String[] args) {
		new RegenieBedAnnot().instanceMainWithExit(args);
	}

}
