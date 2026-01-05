/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.regenie;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLine;
import com.github.lindenb.jvarkit.bed.BedLineCodec;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
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
modificationDate="202050515",
jvarkit_amalgamion = true,
generate_doc = true
)
public class RegenieBedAnnot extends AbstractRegenieAnnot {
	private static final Logger LOG = Logger.of(RegenieBedAnnot.class);
	@Parameter(names = {"-B","--bed"}, description = "custom bed file chrom/start/end/annotation/name[/score].",required = true)
	private Path userBed;
	@Parameter(names = {"-m","--min-length"}, description = "slop each BED records in 5' and 3' so the minimal LENGTH is 'm'. Multiple are comma separated ")
	private String min_len_str="0";
	@Parameter(names = {"--chrom"}, description = "process only that chromosome")
	private String only_that_chrom=null;
	@Parameter(names = {"--noXY"}, description = "skip X/Y chromosome")
	private boolean skipXY=false;	
	private final IntervalTreeMap<List<UserBed>> interval2userbed = new IntervalTreeMap<>();
	
	private int[] min_length_sorted_array = null;
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	private static class UserBed extends SimpleInterval {
		String annotation;
		String gene_name;
		double score;
		UserBed(String ctg,int start,int end) {
			super(ctg,start,end);
		}
		@Override
		public String toString() {
			return super.toString()+" annotation:"+annotation+" gene_name:"+gene_name+" score:"+score;
			}
		}
	
	private static boolean isXY(final String ctg) {
		if(ctg.equals("X") || ctg.equals("Y")) return true;
		if(ctg.equals("chrX") || ctg.equals("chrY")) return true;
		if(ctg.equals("23") || ctg.equals("24")) return true;
		if(ctg.equals("chr23") || ctg.equals("chr24")) return true;
		return false;	
		}
	
	@Override
	protected VCFHeader initVcfHeader(final VCFHeader h) {
		final String default_annotation_value=IOUtils.getFilenameWithoutCommonSuffixes(this.userBed).replaceAll("[^a-zA-Z0-9]+","_");
			
		
		this.min_length_sorted_array = Arrays.stream(CharSplitter.COMMA.split(this.min_len_str)).
				filter(S->!StringUtils.isBlank(S)).
				mapToInt(S->Integer.parseInt(S)).
				filter(S->S>=0).
				sorted().
				toArray();
		if(this.min_length_sorted_array.length==0)  {
			this.min_length_sorted_array=new int[] {0};
			}
		
		final Set<String> allowed_chroms = new HashSet<>();
		if(!StringUtils.isBlank(this.only_that_chrom)) {
			final String ctg=this.only_that_chrom;
			if(!(isXY(ctg) && this.skipXY)) {
				allowed_chroms.add(ctg);
				if(ctg.startsWith("chr")) {
					allowed_chroms.add(ctg.substring(3));
					}
				else
					{
					allowed_chroms.add("chr"+ctg);
					}
				// in REGENIE X and Y are merged
				if(isXY(ctg)) {
					allowed_chroms.add("chrX");
					allowed_chroms.add("X");
					allowed_chroms.add("23");
					allowed_chroms.add("chrY");
					allowed_chroms.add("Y");
					allowed_chroms.add("24");
					}
				}
			}
		final Pattern annotation_regex = Pattern.compile("[A-Za-z][A-Za-z_0-9]*");
		try(BufferedReader br=IOUtils.openPathForBufferedReading(this.userBed)) {
			final BedLineCodec bc = new  BedLineCodec();
			final Map<String,List<UserBed>> gene2intervals=new HashMap<>();
			String line;
			while((line=br.readLine())!=null) {
				if(BedLine.isBedHeader(line)) continue;
				final BedLine rec = bc.decode(line);
				if(rec==null) continue;
				final String ctg = fixContig(rec.getContig());
				if(!allowed_chroms.isEmpty() && !allowed_chroms.contains(ctg)) continue;
				if(this.skipXY && isXY(ctg)) continue;
				
				final UserBed ub = new UserBed(ctg, rec.getStart(), rec.getEnd());
				ub.annotation = StringUtils.ifBlank(rec.get(3),default_annotation_value);
				if(!annotation_regex.matcher(ub.annotation).matches()) {
					throw new IOException("bad annotation in "+line+" should match "+annotation_regex.pattern());
					}
				ub.gene_name = rec.get(4);
				if(StringUtils.isBlank(ub.gene_name) || ub.gene_name.equals(".")) {
					throw new IOException("empty title in bed line "+line);
					}
				
				
				final String scoreStr=rec.getOrDefault(5,"");
				if(StringUtils.isBlank(scoreStr) || scoreStr.equals(".")) {
					ub.score = 1.0;
					}
				else
					{
					ub.score = Double.parseDouble(scoreStr);
					}
				if(ub.score<=0) {
					if(LOG.isDebug()) {
						LOG.debug("skipping "+ub+" because score <=0");
						}
					continue;
					}
				
				List<UserBed> L = gene2intervals.get(ub.gene_name);
				if(L==null) {
					L = new ArrayList<>();
					gene2intervals.put(ub.gene_name, L);
					}
				else
					{
					/* same gene but different annotation at the same place , with same score */
					final Optional<UserBed> opt=L.
						stream().
						filter(R->R.overlaps(ub)).
						filter(UB->!UB.annotation.equals(ub.annotation)).
						filter(UB->UB.score==ub.score).
						findFirst();
					
					if(opt.isPresent()) {
						throw new IOException("for the same gene ("+ ub.gene_name+"), we have two different overlapping regions, with the same score but with different annotations "+ ub+" "+opt.get());
						}
					}
				
				L.add(ub);
				}
			
			
			for(String geneName: gene2intervals.keySet()) {
				final List<UserBed> user_bed_list = gene2intervals.get(geneName);
				List<Locatable> previous_lengths=null;
				for(int min_len : this.min_length_sorted_array) {
					// check extendingh is required or useless by comparing with previous start/end
					if(previous_lengths!=null) {
						boolean same=true;
						for(int i=0;i< user_bed_list.size();++i) {
							final UserBed ub = user_bed_list.get(i);
							final Locatable prev  = previous_lengths.get(i);
							if(ub.getStart()!=prev.getStart() || ub.getEnd()!=prev.getEnd()) {
								same=false;
								break;
								}
							}
						if(same) continue;
						}
					previous_lengths = new ArrayList<>(user_bed_list.size());
					for(int i=0;i< user_bed_list.size();++i) {
						final UserBed ub_src = user_bed_list.get(i);
						int chromStart = ub_src.getStart();
						int chromEnd = ub_src.getEnd();
						if(ub_src.getLengthOnReference()< min_len) {
							final int x2 = Math.max((min_len - ub_src.getLengthOnReference())/2,1);
							chromStart = Math.max(chromStart-x2, 1);
							chromEnd += x2;
							}
						
						//we need to make a copy because we rename and extend user bed
						final UserBed ub = new UserBed(ub_src.getContig(),chromStart,chromEnd);
						ub.gene_name = geneName+(min_len>0?"_x"+min_len:"");
						ub.annotation = ub_src.annotation;
						ub.score = ub_src.score;

						
						final Interval reg  = new Interval(ub);
						previous_lengths.add( reg);
						
						
						List<UserBed> beds = 	this.interval2userbed.get(reg);
						if(beds==null) {
							beds = new ArrayList<>();
							this.interval2userbed.put(reg, beds);
							}
						
							
						beds.add(ub);
						}
					}
				}
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		return super.initVcfHeader(h);
		}
		
	@Override
	protected void dump(final PrintWriter w,final VariantContext ctx) throws Exception {	
		
		final Map<String,List<UserBed>> gene_name_to_user_beds = this.interval2userbed.getOverlapping(ctx).
				stream().
				flatMap(T->T.stream()).
				collect(Collectors.groupingBy(UB->UB.gene_name))
				;
		if(gene_name_to_user_beds.isEmpty()) return;
		
		
		for(String gene_name:gene_name_to_user_beds.keySet()) {
			// sort on best score for this gene
			final UserBed ub = gene_name_to_user_beds.get(gene_name).
					stream().
					sorted((A,B)->Double.compare(B.score,A.score)). //highest score first
					findFirst().
					get();
			
			final Variation v = new Variation();
			v.contig = fixContig(ctx.getContig());
			v.pos = ctx.getStart();
			v.id = makeID(ctx);
			v.gene = ub.gene_name;
			v.prediction = ub.annotation;
			v.score = OptionalDouble.of(ub.score);
			v.cadd = getCaddScore(ctx);
			v.is_singleton = isSingleton(ctx);
			v.frequency = getFrequency(ctx);
			print(w,v);
			}
		}
	

	public static void main(final String[] args) {
		new RegenieBedAnnot().instanceMainWithExit(args);
	}

}
