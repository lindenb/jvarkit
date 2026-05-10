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
package com.github.lindenb.jvarkit.tools.vcf4igv;

import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bio.AcidNucleics;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.dict.SequenceDictionaryExtractor;
import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.SampleSheet;
import com.github.lindenb.jvarkit.io.SampleSheetFactory;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.google.gson.stream.JsonWriter;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;

/**
BEGIN_DOC


## Example
 
 ```
 java -jar dist/jvarkit.jar vcf4igv -R src/test/resources/rotavirus_rf.fa --samplesheet jeter.csv  src/test/resources/rotavirus_rf.vcf.gz | python3 -m json.tool

[
    {
        "fasta": "/home/lindenb/src/jvarkit-git/src/test/resources/rotavirus_rf.fa",
        "interval": "RF01:970-970",
        "contig": "RF01",
        "chromosome": "RF01",
        "start": 970,
        "end": 970,
        "length": 1,
        "ref": "A",
        "alt": "C",
        "description": "RF01:970 CASE/HET:N=0.  CASE/HOM_VAR:N=1 displayed here: S5.  CASE/HOM_REF:N=2 displayed here: S1.  CTRL/HET:N=0.  CTRL/HOM_VAR=0.  CTRL/HOM_REF:N=2 displayed here: S2. ",
        "bams": [
            {
                "sample": "S5",
                "bam": "src/test/resources/S5.bam",
                "status": "case"
            },
            {
                "sample": "S1",
                "bam": "src/test/resources/S1.bam",
                "status": "case"
            },
            {
                "sample": "S2",
                "bam": "src/test/resources/S2.bam",
                "status": "control"
            }
        ]
    },
(...)
```

 END_DOC
 
 */

@Program(
		name="vcf4igv",
		description="Prepare IGV sessions file and json for nextflow/igv_report",
		keywords={"vcf","igv","json"},
		creationDate="20260506",
		modificationDate="20260506",
		jvarkit_amalgamion =  true,
		menu="VCF Manipulation"
		)
public class VcfForIGV extends Launcher {
	private static final Logger LOG = Logger.of( VcfForIGV.class);

	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required = true)
	private Path fasta = null;
	@Parameter(names={"--samplesheet"},description="samplesheet. TSV or CSV file with the following header: 'bam,sample,status'. 'bam' is required, 'status' must be 'case' or 'control'",required=true)
	private Path samplesheet = null;
	@Parameter(names={"--session-dir"},description="if defined, save the IGV session XML files in that directory")
	private Path igv_session_dir = null;
	@Parameter(names={"--count"},description=
			"For Case-controls (at least one 'case' and one 'control' in the samplesheet' ) 8 comma integers representing the maximum number of bam to display: case-HOM_REF,case-HET,case-HOM_VAR,case-NO_CALL,ctrl-HOM_REF,ctrl-HET,ctrl-HOM_VAR,ctr-NO_CALL ."
			+ " A negative number is 'infinite'."
			+ " If there is no case  4 comma separated integers:HOM_REF,HET,HOM_VAR,NO_CALL (other fields are ignored)")
	private String countCaseCtr="1,5,1,0,1,5,1,0";


	public VcfForIGV() {
		
	}
	
	private static class BamInfo {
		String sample;
		String status;
		Path bam;
		SAMFileHeader getHeader() throws IOException {
			final SamReaderFactory srf = SamReaderFactory.make();
			try(SamReader sr=srf.open(this.bam)) {
				return sr.getFileHeader();
				}
			}
		boolean isCase() {
			return this.status.equals("case");
			}
		boolean isControl() {
			return !isCase();
			}
		}
	
	
	static private int compareGT(final Genotype g1, final Genotype g2) {
		int i  = Integer.compare((g1.isFiltered()?1:-1),(g2.isFiltered()?1:-1));
		if(i!=0) return i;
		if(g1.hasGQ() && g2.hasGQ()) {
			i=  Integer.compare(g2.getGQ(),g1.getGQ());
			if(i!=0) return i;
			}
		if(g1.hasDP() && g2.hasDP()) {
			i=  Integer.compare(g2.getDP(),g1.getDP());
			if(i!=0) return i;
			}
		if(g1.hasAD() && g2.hasAD()) {
			final int[] a1=g1.getAD();
			final int[] a2=g2.getAD();
			if(a1.length==2 && a2.length==2) {
				final int n1= a1[0]+a1[1];
				final int n2= a2[0]+a2[1];
				if(n1>0 && n2>0) {
					float f1 = a1[1]/(float)n1;
					float f2 = a2[1]/(float)n2;
					if(g1.isHet() && g2.isHet()) {
						i= Float.compare(Math.abs(0.5f-f1), Math.abs(0.5f-f2));
						if(i!=0) return i;
						}
					else if(g1.isHomRef() && g2.isHomRef()) {
						i= Float.compare(f1,f2);//lowest is bestter
						if(i!=0) return i;
						}
					else if(g1.isHomVar() && g2.isHomVar()) {
						i= Float.compare(f2,f1);//highest is bestter
						if(i!=0) return i;
						}
					}
				}
			}
		return g1.getSampleName().compareTo(g2.getSampleName());
		}
	
	private static List<Genotype> makeList(
			final Map<String,BamInfo> sample2baminfo,
			final VariantContext ctx,
			final boolean is_case,
			final GenotypeType gtype
			) {
			return ctx.getGenotypes()
				.stream()
				.filter(G->sample2baminfo.containsKey(G.getSampleName()) && sample2baminfo.get(G.getSampleName()).isCase()== is_case)
				.filter(G->G.getType().equals(gtype))
				.sorted(VcfForIGV::compareGT)
				.collect(Collectors.toCollection(ArrayList::new));
			}
	
	private static List<Genotype> makeList(
			final Map<String,BamInfo> sample2baminfo,
			final VariantContext ctx,
			final GenotypeType gtype
			) {
			return ctx.getGenotypes()
				.stream()
				.filter(G->sample2baminfo.containsKey(G.getSampleName()))
				.filter(G->G.getType().equals(gtype))
				.sorted(VcfForIGV::compareGT)
				.collect(Collectors.toCollection(ArrayList::new));
			}
	
	private static List<Genotype> limit(final List<Genotype> L,int N) {
		if(N<0) return L;
		while(L.size()>N) L.remove(L.size()-1);
		return L;
		}
	private static String toAbsolutePath(Path p) {
		try {
			final Path p2 = p.toRealPath();//resolve symlink for nextflow
			p=p2;
			}
		catch(Throwable err) {
			
			}
		
		try {
			return p.toAbsolutePath().toString();
		}
		catch(Throwable err) {
			return p.toString();
		}
	}
	
	@Override
	public int doWork(final List<String> args) {
		final int COL_HOM_REF=0;
		final int COL_HET=1;
		final int COL_HOM_VAR=2;
		final int COL_NO_CALL=3;
		
		try {
			final SAMSequenceDictionary dict = new SequenceDictionaryExtractor().extractRequiredDictionary(fasta);
			
			final Optional<String> buildName= SequenceDictionaryUtils.getBuildName(dict);
			
			final SampleSheet sampleSheet =  new SampleSheetFactory().
					splitter(CharSplitter.forFilename(this.samplesheet.getFileName().toString())).
					of(this.samplesheet);
			sampleSheet.getHeader().assertColumnExists("bam");
			final Map<String,BamInfo> sample2baminfo = new HashMap<>();
			
			final XMLOutputFactory xmlOutputFactory = XMLOutputFactory.newFactory();
			
			
			final String input = super.oneFileOrNull(args);
			try(VCFIterator r = input==null?new VCFIteratorBuilder().open(System.in):new VCFIteratorBuilder().open(input)) {
				final VCFHeader h = r.getHeader();
				final SAMSequenceDictionary vcfdict = h.getSequenceDictionary();
				if(vcfdict!=null ) {
					SequenceUtil.assertSequenceDictionariesEqual(dict, vcfdict);
					}
				final Set<String> samples_in_vcf= h.getSampleNameToOffset().keySet();
				
				for(FileHeader.RowMap rowMap : sampleSheet ) {
					final BamInfo bi = new BamInfo();
					if(StringUtils.isBlank(rowMap.get("bam"))) {
						LOG.info("skipping ("+rowMap+") because bam is empty");
						continue;
						}
					bi.bam = Paths.get(rowMap.get("bam"));
					if(!Files.exists(bi.bam)) {
						LOG.info("skipping ("+rowMap+") because cannot find bam \""+bi.bam+"\"");
						continue;
						}
					final SAMFileHeader samH  = bi.getHeader();
					final SAMSequenceDictionary samdict =  SequenceDictionaryUtils.extractRequired(samH);
					if(!SequenceUtil.areSequenceDictionariesEqual(dict, samdict)) {
						LOG.info("skipping ("+rowMap+") is not mapped of "+this.fasta);
						continue;
						}
						
					bi.sample= rowMap.getOrDefault("sample","");
					if(StringUtils.isBlank(bi.sample)) {
						bi.sample = samH.getReadGroups().stream().map(RG->RG.getSample()).filter(S->!StringUtils.isBlank(S)).findFirst().orElse("");
						}
					
					if(!samples_in_vcf.contains(bi.sample)) {
						LOG.info("skipping ("+rowMap+") because cannot find sample "+bi.sample+" in vcf");
						continue;
						}
					if(sample2baminfo.containsKey(bi.sample)) {
						LOG.info("skipping "+rowMap+" because it was defined twice.");
						continue;
						}
					bi.status = StringUtils.ifBlank( rowMap.getOrDefault("status", "control"),"control");
					if(!(bi.status.equals("case") || bi.status.equals("control"))) {
						LOG.error("status must be case or control (default) in "+rowMap);
						return -1;
						}
					
					
					sample2baminfo.put(bi.sample, bi);
					}
				
				if(sample2baminfo.isEmpty()) {
					LOG.error("no samples in common between vcf and "+this.samplesheet);
					return -1;
					}
				final boolean is_case_controls = 
						sample2baminfo.values().stream().anyMatch(BI->BI.isCase()) &&
						sample2baminfo.values().stream().anyMatch(BI->BI.isControl())
						;
				final int counts[]=  Arrays.stream(CharSplitter.COMMA.split(this.countCaseCtr))
						.mapToInt(S->Integer.parseInt(S))
						.toArray();
				if(is_case_controls) {
					if(counts.length!=8) throw new JvarkitException.TokenErrors(8, CharSplitter.COMMA.split(this.countCaseCtr));
					}
				else
					{
					// just use 3 first
					if(!(counts.length!=4 || counts.length!=8)) throw new JvarkitException.TokenErrors(4, CharSplitter.COMMA.split(this.countCaseCtr));
					}
				
				try(JsonWriter jw = new JsonWriter(super.openPathOrStdoutAsPrintWriter(this.outputFile)) ) {
					jw.beginArray();
					
					while(r.hasNext()) {
						final VariantContext ctx = r.next();
						if(dict.getSequence(ctx.getContig())==null) {
							LOG.warning("skpipping "+ctx.getContig()+":"+ctx.getStart()+" because contig is not in fasta ref");
							continue;
							}
						final StringBuilder description = new StringBuilder();
						description.append(ctx.getContig()+":"+ctx.getStart());
						if(ctx.getLengthOnReference()!=1) {
							description.append("-"+ctx.getEnd());
							description.append(" length:"+StringUtils.niceInt(ctx.getLengthOnReference()));
							}
						
						if(buildName.isPresent()) {
							description.append(" ").append(buildName.get());
							}
						if(ctx.hasAttribute(VCFConstants.SVTYPE)) {
							description.append(" SVTYPE:").append(ctx.getAttributeAsString(VCFConstants.SVTYPE,"."));
							}
						if(ctx.hasAttribute("SVLEN")) {
							description.append(" SVLEN:").append(ctx.getAttribute("SVLEN"));
							}
						
						final List<Genotype> L=new ArrayList<>();
						if(is_case_controls) {
							final List<Genotype> case_het =makeList(sample2baminfo,ctx,true,GenotypeType.HET);
							description.append(" CASE/HET:N=").append(case_het.size());
							limit(case_het,counts[COL_HET]);
							if(!case_het.isEmpty()) {
								description
									.append(" displayed here: ")
									.append(case_het.stream().map(G->G.getSampleName()).collect(Collectors.joining(",")));
								}
							description.append(". ");
							
							final List<Genotype> case_hom_var =makeList(sample2baminfo,ctx,true,GenotypeType.HOM_VAR);
							description.append(" CASE/HOM_VAR:N=").append(case_hom_var.size());
							limit(case_hom_var,counts[COL_HOM_VAR]);
							if(!case_hom_var.isEmpty()) {
								description
									
									.append(" displayed here: ")
									.append(case_hom_var.stream().map(G->G.getSampleName()).collect(Collectors.joining(",")));
								}
							description.append(". ");
							
							final List<Genotype> case_homref =makeList(sample2baminfo,ctx,true,GenotypeType.HOM_REF);
							description.append(" CASE/HOM_REF:N=").append(case_homref.size());
							limit(case_homref,counts[COL_HOM_REF]);
							if(!case_homref.isEmpty()) {
								description
									
									.append(" displayed here: ")
									.append(case_homref.stream().map(G->G.getSampleName()).collect(Collectors.joining(",")));
								}
							description.append(". ");
							
							
							final List<Genotype> case_nocall =makeList(sample2baminfo,ctx,true,GenotypeType.NO_CALL);
							description.append(" CASE/NO_CALL:N=").append(case_nocall.size());
							limit(case_nocall,counts[COL_NO_CALL]);
							if(!case_nocall.isEmpty()) {
								description
									.append(" displayed here: ")
									.append(case_nocall.stream().map(G->G.getSampleName()).collect(Collectors.joining(",")));
								}
							description.append(". ");
							
							
							final List<Genotype> ctrl_het =makeList(sample2baminfo,ctx,false,GenotypeType.HET);
							description.append(" CTRL/HET:N=").append(ctrl_het.size());
							limit(ctrl_het,counts[COL_HET+4]);
							if(!ctrl_het.isEmpty()) {
								description
									.append(" displayed here: ")
									.append(ctrl_het.stream().map(G->G.getSampleName()).collect(Collectors.joining(",")));
								}
							description.append(". ");
							
							final List<Genotype> ctrl_hom_var =makeList(sample2baminfo,ctx,false,GenotypeType.HOM_VAR);
							description.append(" CTRL/HOM_VAR=").append(ctrl_hom_var.size());
							limit(ctrl_hom_var,counts[COL_HOM_VAR+4]);
							if(!ctrl_hom_var.isEmpty()) {
								description
									.append(" displayed here: ")
									.append(ctrl_hom_var.stream().map(G->G.getSampleName()).collect(Collectors.joining(",")));
								}
							description.append(". ");
							
							final List<Genotype> ctrl_homref =makeList(sample2baminfo,ctx,false,GenotypeType.HOM_REF);
							description.append(" CTRL/HOM_REF:N=").append(ctrl_homref.size());
							limit(ctrl_homref,counts[COL_HOM_REF+4]);
							if(!ctrl_homref.isEmpty()) {
								description
									.append(" displayed here: ")
									.append(ctrl_homref.stream().map(G->G.getSampleName()).collect(Collectors.joining(",")));
								}
							description.append(". ");
							
							final List<Genotype> ctrl_nocall =makeList(sample2baminfo,ctx,false,GenotypeType.NO_CALL);
							description.append(" CTRL/NO_CALL:N=").append(ctrl_nocall.size());
							limit(ctrl_nocall,counts[COL_NO_CALL+4]);
							if(!ctrl_nocall.isEmpty()) {
								description
									.append(" displayed here: ")
									.append(ctrl_nocall.stream().map(G->G.getSampleName()).collect(Collectors.joining(",")));
								}
							description.append(". ");
							
							
							L.addAll(case_het);
							L.addAll(case_hom_var);
							L.addAll(case_homref);
							L.addAll(case_nocall);
							L.addAll(ctrl_het);
							L.addAll(ctrl_hom_var);
							L.addAll(ctrl_homref);
							L.addAll(ctrl_nocall);
							}
						else
							{
							final List<Genotype> any_het = makeList(sample2baminfo,ctx,GenotypeType.HET);
							description.append(" HET:N=").append(any_het.size());
							limit(any_het,counts[COL_HET]);
							if(!any_het.isEmpty()) {
								description
									.append(" displayed here: ")
									.append(any_het.stream().map(G->G.getSampleName()).collect(Collectors.joining(",")));
								}
							description.append(". ");
							
							final List<Genotype> any_hom_var = makeList(sample2baminfo,ctx,GenotypeType.HOM_VAR);
							description.append(" HOM_VAR:N=").append(any_hom_var.size());
							limit(any_hom_var,counts[COL_HOM_VAR]);
							if(!any_hom_var.isEmpty()) {
								description
									.append(" displayed here: ")
									.append(any_hom_var.stream().map(G->G.getSampleName()).collect(Collectors.joining(",")));
								}
							description.append(". ");
							
							final List<Genotype> any_homref = makeList(sample2baminfo,ctx,GenotypeType.HOM_REF);
							description.append(" HOM_REF:N=").append(any_homref.size());
							limit(any_homref,counts[COL_HOM_REF]);
							if(!any_homref.isEmpty()) {
								description
									.append(" displayed here: ")
									.append(any_homref.stream().map(G->G.getSampleName()).collect(Collectors.joining(",")));
								}
							description.append(". ");
							
							final List<Genotype> any_nocall = makeList(sample2baminfo,ctx,GenotypeType.NO_CALL);
							description.append(" NO_CALL:N=").append(any_nocall.size());
							limit(any_nocall,counts[COL_NO_CALL]);
							if(!any_nocall.isEmpty()) {
								description
									.append(" displayed here: ")
									.append(any_nocall.stream().map(G->G.getSampleName()).collect(Collectors.joining(",")));
								}
							description.append(". ");
							
							L.addAll(any_homref);
							L.addAll(any_het);
							L.addAll(any_hom_var);
							L.addAll(any_nocall);
							}
						if(L.isEmpty()) continue;
						
						jw.beginObject();
						jw.name("fasta");jw.value(toAbsolutePath(fasta));
						jw.name("interval");jw.value(ctx.getContig()+":"+ctx.getStart()+"-"+ctx.getEnd());
						jw.name("contig");jw.value(ctx.getContig());
						jw.name("chromosome");jw.value(ctx.getContig());
						jw.name("start");jw.value(ctx.getStart());
						jw.name("end");jw.value(ctx.getEnd());
						jw.name("length");jw.value(ctx.getLengthOnReference());

						if(ctx.hasID()) {
							jw.name("variant-id");jw.value(ctx.getID());
							}
						
						jw.name("filtered");jw.value(ctx.isFiltered());
						
						if(AcidNucleics.isATGC(ctx.getReference())) {
							jw.name("ref");jw.value(ctx.getReference().getDisplayString());
							}
						if(ctx.getNAlleles()==2 && AcidNucleics.isATGC(ctx.getAlleles().get(1))) {
							jw.name("alt");jw.value(ctx.getAlleles().get(1).getDisplayString());
							}
						if(ctx.hasAttribute(VCFConstants.SVTYPE)) {
							jw.name("svtype");
							jw.value(ctx.getAttributeAsString(VCFConstants.SVTYPE,"."));
							}
						jw.name("description");jw.value(description.toString());

						jw.name("bams");
						jw.beginArray();
						for(Genotype gt:L) {
							final BamInfo bi = sample2baminfo.get(gt.getSampleName());
							jw.beginObject();
							jw.name("sample");jw.value(bi.sample);
							jw.name("bam");jw.value(bi.bam.toString());
							jw.name("status");jw.value(bi.status);
							jw.name("description");jw.value(gt.getType().name()
									+(gt.hasGQ()?" GQ:"+gt.getGQ():"")
									+(gt.hasDP()?" DP:"+gt.getDP():"")
									+(gt.hasAD()?" AD:"+Arrays.stream(gt.getAD()).mapToObj(D->String.valueOf(D)).collect(Collectors.joining(",")):"")
									);
							jw.endObject();
							}
						
						jw.endArray();
						
						jw.endObject();
						
						if(this.igv_session_dir!=null) {
							final String fname = IOUtils.escapePath(ctx.getContig())+"_"+ ctx.getStart()+"_"+ctx.getEnd()+"_"+StringUtils.md5(ctx.getAlleles().stream().map(S->S.toString()).collect(Collectors.joining(",")))+".igv_session.xml";
							final Path session_path = this.igv_session_dir.resolve(fname);
							try(OutputStream os = Files.newOutputStream(session_path)) {
								final XMLStreamWriter w = xmlOutputFactory.createXMLStreamWriter(os, "UTF-8");
								w.writeStartDocument("UTF-8", "1.0");
								w.writeStartElement("Session");
								w.writeAttribute("genome",toAbsolutePath(this.fasta));
								w.writeAttribute("locus",ctx.getContig()+":"+ctx.getStart()+"-"+ctx.getEnd());
								w.writeAttribute("version","8");
								w.writeStartElement("Resources");
								for(BamInfo bi: L.stream().map(G->sample2baminfo.get(G.getSampleName())).collect(Collectors.toList())) {
									w.writeEmptyElement("Resource");
									w.writeAttribute("path",toAbsolutePath(bi.bam));
									if(bi.bam.getFileName().toString().endsWith(FileExtensions.CRAM)) {
										w.writeAttribute("type","cram");
										}
									else if(bi.bam.getFileName().toString().endsWith(FileExtensions.BAM)) {
										w.writeAttribute("type","bam");
										}
									}
								w.writeEndElement();//end Resources

								w.writeEndElement();//end Session
								w.writeEndDocument();
								w.flush();
								os.flush();
								}
							}
						}// while r.hasNext
					
					jw.endArray();
					jw.flush();
					}
				}//try
			return 0;
		    }
		catch(final Throwable err ) {
			LOG.error(err);
			return -1;
		}
	}

	public static void main(String[] args) {
		new VcfForIGV().instanceMainWithExit(args);
	}
}
