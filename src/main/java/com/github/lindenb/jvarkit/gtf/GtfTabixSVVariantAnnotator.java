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
package com.github.lindenb.jvarkit.gtf;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.github.lindenb.jvarkit.lang.AttributeMap;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.tabix.AbstractTabixVariantAnnotator;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFCodec;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFLine;

import htsjdk.samtools.util.CoordMath;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/*** tabix annotator using GTF to annotate SV */
public class GtfTabixSVVariantAnnotator extends AbstractTabixVariantAnnotator {
	public static final String OPT_DESC="GTF file";
	public static final String EXTEND_KEY ="extend";
	
	private int extend;
	private VCFInfoHeaderLine hdrInCDS;
	private VCFInfoHeaderLine hdrInExon;
	private VCFInfoHeaderLine hdrInTranscript;
	private VCFInfoHeaderLine hdrInGene;
	private VCFInfoHeaderLine hdrInBioTypeProteinCoding;
	
	private VCFInfoHeaderLine hdrGeneIDSet;
	private VCFInfoHeaderLine hdrGeneNameSet;
	private VCFInfoHeaderLine hdrGeneBiotypeSet;

	private VCFInfoHeaderLine hdrIntronDelTranscriptIds;
	private VCFInfoHeaderLine hdrIntronDelFlag;

	
	private VCFInfoHeaderLine hdrInDownstream;
	private VCFInfoHeaderLine hdrInUpstream;
	private final GTFCodec gtfCodec;
	public GtfTabixSVVariantAnnotator(final String uri, AttributeMap meta) throws IOException {
		super(uri);
		if(meta==null) meta = AttributeMap.empty();
		this.extend = meta.getIntAttribute(EXTEND_KEY).orElse(0);
		this.gtfCodec = new GTFCodec();
		}

	@Override
	public void fillHeader(final VCFHeader header) {
		if(super.tabixReader==null) return;
		String helpCmd = " in " + getUri();
		this.hdrInCDS = new VCFInfoHeaderLine("IN_CDS",1, VCFHeaderLineType.Flag, "SV was found overlaping a CDS  "+ helpCmd );
		header.addMetaDataLine(this.hdrInCDS);
		this.hdrInExon = new VCFInfoHeaderLine("IN_EXON",1, VCFHeaderLineType.Flag, "SV was found overlaping an exon "+ helpCmd );
		header.addMetaDataLine(this.hdrInExon);
		this.hdrInTranscript = new VCFInfoHeaderLine("IN_TRANSCRIPT",1, VCFHeaderLineType.Flag, "SV was found overlaping a transcript "+helpCmd);
		header.addMetaDataLine(this.hdrInTranscript);
		this.hdrInGene = new VCFInfoHeaderLine("IN_GENE",1, VCFHeaderLineType.Flag, "SV was found overlaping a gene "+helpCmd);
		header.addMetaDataLine(this.hdrInGene);
		this.hdrInBioTypeProteinCoding = new VCFInfoHeaderLine("IN_PROTEIN_CODING",1, VCFHeaderLineType.Flag, "SV was found overlaping a gene with (gene_biotype|transcript_type)=protein_coding "+helpCmd);
		header.addMetaDataLine(this.hdrInBioTypeProteinCoding);

		
		this.hdrGeneIDSet = new VCFInfoHeaderLine("GENE_ID",VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene IDs "+helpCmd);
		header.addMetaDataLine(this.hdrGeneIDSet);
		this.hdrGeneNameSet = new VCFInfoHeaderLine("GENE_NAME",VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene Names "+helpCmd);
		header.addMetaDataLine(this.hdrGeneNameSet);
		this.hdrGeneBiotypeSet = new VCFInfoHeaderLine("GENE_BIOTYPE",VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Gene Biotypes "+helpCmd);
		header.addMetaDataLine(this.hdrGeneBiotypeSet);

		
		// intron skip
		this.hdrIntronDelFlag = new VCFInfoHeaderLine("INTRON_SKIP",1, VCFHeaderLineType.Flag, "SV has the coordinate of an intron (retrogene...) "+helpCmd);
		header.addMetaDataLine(this.hdrIntronDelFlag);
		this.hdrIntronDelTranscriptIds = new VCFInfoHeaderLine("INTRON_SKIP_TRANSCRIPTS",VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Transcript ids where SV has the coordinate of an intron (retrogene...) "+helpCmd);
		header.addMetaDataLine(this.hdrIntronDelTranscriptIds);

		
		helpCmd+=" within a distance of "+this.extend+" bp.";
		
		this.hdrInDownstream = new VCFInfoHeaderLine("DOWNSTREAM_GENE",VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "SV  found overlaping downstream of a gene  "+helpCmd);
		header.addMetaDataLine(this.hdrInDownstream);
		this.hdrInUpstream = new VCFInfoHeaderLine("UPSTREAM_GENE",VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "SV  found overlaping upstream of a gene "+helpCmd);
		header.addMetaDataLine(this.hdrInUpstream);
		}
	@Override
	public List<VariantContext> annotate(final VariantContext ctx) throws IOException  {
		boolean in_cds_flag = false;
		boolean in_exon_flag = false;
		boolean in_transcript_flag = false;
		boolean in_gene_flag = false;
		boolean in_protein_coding = false;
		if(super.tabixReader ==null || !hasContig(ctx)) return Collections.singletonList(ctx);
		int extend_start = Math.max(0, ctx.getStart()- this.extend);
		int extend_end = ctx.getEnd() + this.extend;
		final Set<String> downstream_genes = new HashSet<>();
		final Set<String> upstream_genes = new HashSet<>();
		
		final Set<String> all_genes_names = new HashSet<>();
		final Set<String> all_genes_ids = new HashSet<>();
		final Set<String> all_genes_biotypes = new HashSet<>();
		TabixReader.Iterator r = this.tabixReader.query(contig(ctx), extend_start, extend_end);
		for(;;) {
			String line = r.next();
			if(line==null) break;
			final String tokens[] = CharSplitter.TAB.split(line);
			final int rec_start = Integer.parseInt(tokens[3]);
			final int rec_end = Integer.parseInt(tokens[4]);
			final GTFLine gtfRecord = this.gtfCodec.decode(tokens);
			if(CoordMath.overlaps(ctx.getStart(), ctx.getEnd(), rec_start, rec_end)) {
				if(gtfRecord.isCDS()) {in_cds_flag=true; in_exon_flag=true; in_transcript_flag=true; in_gene_flag=true;}
				else if(gtfRecord.isExon()) {in_exon_flag=true; in_transcript_flag=true; in_gene_flag=true;}
				else if(gtfRecord.isTranscript()) {in_transcript_flag=true; in_gene_flag=true;}
				else if(gtfRecord.isGene()) {in_gene_flag=true;}
				if(gtfRecord.hasAttribute("gene_id")) {
					all_genes_ids.add(gtfRecord.getAttributes().get("gene_id"));
					}
				if(gtfRecord.hasAttribute("gene_name")) {
					all_genes_names.add(gtfRecord.getAttributes().get("gene_name"));
					}
				if(gtfRecord.hasAttribute("gene_biotype")) {
					all_genes_biotypes.add(gtfRecord.getAttributes().get("gene_biotype"));
					}
				if(gtfRecord.getAttributes().entrySet().stream().
					filter(KV->KV.getValue().equals("protein_coding")).
					map(KV->KV.getKey()).
					anyMatch(KEY->KEY.equals("transcript_type") || KEY.equals("gene_type") || KEY.equals("transcript_biotype") || KEY.equals("gene_biotype"))) {
					in_protein_coding = true;
					}
				}
			
			
			
			if(this.extend>0 &&
				!CoordMath.overlaps(ctx.getStart(), ctx.getEnd(), rec_start, rec_end) && 
				CoordMath.overlaps(ctx.getStart(), ctx.getEnd(), extend_start, extend_end) && 
				tokens[2].equals("gene")) {
					final boolean is_plus_strand = gtfRecord.isPostiveStrand();
					boolean is_upstream;
					int distance;
					if(rec_end < ctx.getStart()) {
						is_upstream = is_plus_strand;
						distance = CoordMath.getLength(rec_end, ctx.getStart());
						}
					else if(rec_start > ctx.getEnd()) {
						is_upstream = !is_plus_strand;
						distance = CoordMath.getLength(ctx.getEnd(),rec_start);
						}
					else
						{
						throw new IllegalStateException();
						}
					
					(is_upstream?upstream_genes:downstream_genes).add(
						String.join("|",
						(is_upstream?"upstream":"downstream"),
						gtfRecord.getAttributes().getOrDefault("gene_id","."),
						gtfRecord.getAttributes().getOrDefault("gene_name","."),
						gtfRecord.getAttributes().getOrDefault("gene_biotype","."),
						String.valueOf(distance)
						));
					}
			}
		//seach for intron deletion
		final Set<String> intron_skip_transcript_ids;
		if(in_gene_flag && ctx.getAttribute(VCFConstants.SVTYPE,"").equals("DEL")) {
			final Set<String> transcript_5 = new HashSet<>();
			final Set<String> transcript_3 = new HashSet<>();
			for(int side=0;side<2;side++) {
				final int junction_extend=5;
				final int junction_pos = (side==0?ctx.getStart():ctx.getEnd());
				final int junction_start = Math.max(1, junction_pos-junction_extend);
				final int junction_end = Math.max(1, junction_pos-junction_extend);
				r = this.tabixReader.query(contig(ctx), junction_start ,junction_end);
				for(;;) {
					final String line = r.next();
					if(line==null) break;
					final String tokens[] = CharSplitter.TAB.split(line);
					if(!tokens[2].equals("exon")) continue;
					final GTFLine gtfRecord = this.gtfCodec.decode(tokens);
					if(!CoordMath.overlaps(junction_start, junction_end, gtfRecord.getStart(), gtfRecord.getEnd())) continue;
					final String transcript_id = gtfRecord.getAttributes().getOrDefault("transcript_id","");
					if(StringUtils.isBlank(transcript_id)) continue;
					(side==0?transcript_5:transcript_3).add(transcript_id);
 					}
				}
			transcript_5.retainAll(transcript_3);
			intron_skip_transcript_ids = transcript_5;
			}
		else
			{
			intron_skip_transcript_ids = Collections.emptySet();
			}
		
		final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
		if(in_cds_flag)  vcb.attribute(this.hdrInCDS.getID(), true);
		if(in_exon_flag)  vcb.attribute(this.hdrInExon.getID(), true);
		if(in_transcript_flag)  vcb.attribute(this.hdrInTranscript.getID(), true);
		if(in_gene_flag)  vcb.attribute(this.hdrInGene.getID(), true);
		if(in_protein_coding)  vcb.attribute(this.hdrInBioTypeProteinCoding.getID(), true);
		if(!upstream_genes.isEmpty()) vcb.attribute(this.hdrInUpstream.getID(), new ArrayList<>(upstream_genes));
		if(!downstream_genes.isEmpty()) vcb.attribute(this.hdrInDownstream.getID(), new ArrayList<>(downstream_genes));
		
		if(!all_genes_names.isEmpty()) {
			 vcb.attribute(this.hdrGeneNameSet.getID(), new ArrayList<>(all_genes_names));
			}
		if(!all_genes_biotypes.isEmpty()) {
			 vcb.attribute(this.hdrGeneBiotypeSet.getID(), new ArrayList<>(all_genes_biotypes));
			}
		if(!all_genes_ids.isEmpty()) {
			 vcb.attribute(this.hdrGeneIDSet.getID(), new ArrayList<>(all_genes_ids));
			}
		if(!intron_skip_transcript_ids.isEmpty()) {
			vcb.attribute(this.hdrIntronDelFlag.getID(), true);
			vcb.attribute(this.hdrIntronDelTranscriptIds.getID(), new ArrayList<>(intron_skip_transcript_ids));
			}
		return Collections.singletonList(vcb.make());
		}
	}
