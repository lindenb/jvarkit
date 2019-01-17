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
package com.github.lindenb.jvarkit.tools.lumpysv;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

/**
BEGIN_DOC

## Example

```
 java -jar dist/lumpyvcf2circos.jar  --minsu 50 -inv -bnb -dup  -o tmp  LumpyExpress.vcf.gz \
  && (cd tmp; /path/to/bin/circos  -outputdir . -conf lumpy.circos.conf  )
```

END_DOC
*/
@Program(name="lumpyvcf2circos",
		description="Lumpy to Circos",
		keywords={"lumpy","circos","sv","vcf"}
		)
public class LumpyVcfToCircos extends Launcher {
	
	private static final Logger LOG = Logger.build(LumpyVcfToCircos.class).make();
	
	@Parameter(names={"-o","--output"},description="output directory or zip file")
	private File outputFile = null;
	@Parameter(names={"-p","--prefix"},description="file prefix")
	private String filePrefix = "lumpy.";
	@Parameter(names={"-inv","--inv","--invertions"},description="hide <INV>")
	private boolean hide_inv = false;
	@Parameter(names={"-del","--del","--deletions"},description="hide <DEL>")
	private boolean hide_del = false;
	@Parameter(names={"-dup","--dup","--duplications"},description="hide <DUP>")
	private boolean hide_dup = false;
	@Parameter(names={"-bnb","-bnd","--bnd"},description="hide SVTYPE=BND")
	private boolean hide_bnd = false;
	@Parameter(names={"--su","--minsu"},description="Min supporting reads.")
	private int minSupportingReads = 20;

	
	private static class AbstractSV
		{
		String contig;
		int start;
		int end;
		/** Number of paired-end reads supporting the variant */
		int pe=0;
		/** Number of split-reads supporting the variant */
		int sr=0;
		/** Number of pieces of evidence supporting the variant */
		int su=0;
		void visit(Genotype gt)
			{
			su = gt.getAttributeAsInt("SU",0);	
			sr = gt.getAttributeAsInt("SR",0);	
			pe = gt.getAttributeAsInt("PE",0);	
			}
		void print(PrintWriter pw)
			{
			pw.print(contig);
			pw.print("\t");
			pw.print(start);
			pw.print("\t");
			pw.print(end);
			pw.print("\t");
			pw.print(su);
			pw.println();
			}
		}
	
	private static class SVDel extends AbstractSV
		{
		}
	private static class SVDup extends AbstractSV
		{
		}
	private static class SVInv extends AbstractSV
		{
		}
	private static class SVBnB extends AbstractSV
		{
		String bnb_contig;
		int bnb_start;
		int bnb_end;
		
		@Override void print(PrintWriter pw)
			{
			pw.print(contig);
			pw.print("\t");
			pw.print(start);
			pw.print("\t");
			pw.print(end);
			pw.print("\t");
			pw.print(bnb_contig);
			pw.print("\t");
			pw.print(bnb_start);
			pw.print("\t");
			pw.print(bnb_end);
			pw.println();
			}
		}
	
	private class SampleWriters
		{
		final String sampleName;
		final List<SVDel> deletions=new ArrayList<>();
		final List<SVDup> duplications=new ArrayList<>();
		final List<SVInv> invertions=new ArrayList<>();
		final List<SVBnB> bdns=new ArrayList<>();
		
		SampleWriters(final String sampleName) throws IOException
			{
			this.sampleName = sampleName;
			}
		String createFilename(String s) {
			return LumpyVcfToCircos.this.createFilename(this.sampleName+"."+s);
		}

		String getDeletionFilename() { return this.createFilename("del.txt");}
		String getDuplicationFilename() { return this.createFilename("dup.txt");}
		String getInvertionFilename() { return this.createFilename("inv.txt");}
		String getBnbFilename() { return this.createFilename("bnb.txt");}
		
		public int getNumberOfTracks()
			{
			int n=0;
			n+=(deletions.isEmpty()?0:1);
			n+=(duplications.isEmpty()?0:1);
			n+=(invertions.isEmpty()?0:1);
			n+=(bdns.isEmpty()?0:1);
			return n;
			}
		
		
		
		void visit(final VariantContext vc)
			{
			final List<Allele> alts = vc.getAlternateAlleles();
			if(alts.size()!=1)
				{
				LOG.warning("unknown situtation "+vc);
				return;
				}
			final Genotype genotype = vc.getGenotype(sampleName);
			final Allele alt = alts.get(0);
			if(alt.equals(LumpyConstants.DEL))
				{
				if(LumpyVcfToCircos.this.hide_del) return;
 				final SVDel sv = new SVDel();
				sv.visit(genotype);
				if(sv.su< LumpyVcfToCircos.this.minSupportingReads) return;
				
				sv.contig  = vc.getContig();
				sv.start= vc.getStart();
				sv.end = vc.getAttributeAsInt("END", -1);
				deletions.add(sv);
				}
			else if(alt.equals(LumpyConstants.INV))
				{
				if(LumpyVcfToCircos.this.hide_inv) return;

				final SVInv sv = new SVInv();
				sv.visit(genotype);
				if(sv.su< LumpyVcfToCircos.this.minSupportingReads) return;
				
				sv.contig  = vc.getContig();
				sv.start= vc.getStart();
				sv.end = vc.getAttributeAsInt("END", -1);
				invertions.add(sv);
				}
			else if(alt.equals(LumpyConstants.DUP))
				{
				if(LumpyVcfToCircos.this.hide_dup) return;
				final SVDup sv = new SVDup();
				sv.visit(genotype);
				if(sv.su< LumpyVcfToCircos.this.minSupportingReads) return;
				
				sv.contig  = vc.getContig();
				sv.start= vc.getStart();
				sv.end = vc.getAttributeAsInt("END", -1);
				duplications.add(sv);
				}
			else if(alt.getDisplayString().contains("[") || alt.getDisplayString().contains("]"))
				{
				if(LumpyVcfToCircos.this.hide_bnd) return;
				if(vc.getID().endsWith("_2")) return;
				if(!vc.getID().endsWith("_1"))
					{
					LOG.warning("unknown situtation "+vc);
					return;
					}
				
				
				String display= alt.getDisplayString();
				final SVBnB sv = new SVBnB();
				sv.visit(genotype);
				if(sv.su< LumpyVcfToCircos.this.minSupportingReads) return;
				sv.contig  = vc.getContig();
				sv.start= vc.getStart();
				sv.end = vc.getStart()+1;
				
				int idx0 = display.indexOf("[");
				int idx1 = display.indexOf("]");
				int idx;
				if(idx0==-1) { idx=idx1;}
				else if(idx1==-1) { idx=idx0;}
				else { idx=Math.min(idx0, idx1);}
				
				display = display.substring(idx+1);
				
				idx0 = display.indexOf("[");
				idx1 = display.indexOf("]");
				if(idx0==-1) { idx=idx1;}
				else if(idx1==-1) { idx=idx0;}
				else { idx=Math.min(idx0, idx1);}
				display = display.substring(0,idx);
				
				int colon = display.indexOf(":");
				if(colon==-1)
					{
					LOG.warning("unknown situtation "+vc);
					return;
					}
				
				
				sv.bnb_contig=display.substring(0,colon);
				sv.bnb_start=Integer.parseInt(display.substring(colon+1));
				sv.bnb_end=sv.bnb_start;
				bdns.add(sv);
				}
			else
				{
				LOG.warning("unknown situtation "+vc);
				}
			}
		
		}
	
	private String createFilename(final String s)
		{
		return this.filePrefix+(filePrefix.endsWith(".") && !filePrefix.isEmpty()?"":".")+s;
		}
	
	@Override
	public int doWork(List<String> args) {
	VCFIterator r=null;
	ArchiveFactory archiveFactory=null;
	PrintWriter pw=null;
	PrintWriter conf=null;
	try {
		r = super.openVCFIterator(oneFileOrNull(args));
		
		final VCFHeader header =r.getHeader();
		final SAMSequenceDictionary dict = header.getSequenceDictionary();
		if(dict==null)
			{
			LOG.error("No dictionary in vcf");
			return -1;
			}
		archiveFactory = ArchiveFactory.open(outputFile);
		
		conf=archiveFactory.openWriter(createFilename("ideogram.conf"));
		conf.println("<ideogram>");
		conf.println("show = yes");
		conf.println("<spacing>");
		conf.println("default = 1u");
		conf.println("<pairwise \"/hs/ /hs/\">");
		conf.println("spacing = 0.5r");
		conf.println("</pairwise>");
		conf.println("");
		conf.println("</spacing>");
		conf.println("");
		conf.println("thickness         = 25p");
		conf.println("");
		conf.println("#stroke_thickness = 1");
		conf.println("#stroke_color     = black");
		conf.println("");
		conf.println("fill           = yes");
		conf.println("fill_color     = black");
		conf.println("");
		conf.println("radius         = 1r");
		conf.println("show_label     = yes");
		conf.println("label_font     = default");
		conf.println("label_radius   = dims(ideogram,radius_outer) + 10p");
		conf.println("label_size     = 12p");
		conf.println("label_parallel = yes");
		conf.println("");
		conf.println("show_bands            = yes");
		conf.println("fill_bands            = yes");
		conf.println("band_stroke_thickness = 0");
		conf.println("band_stroke_color     = black");
		conf.println("band_transparency     = 4");
		conf.println("");
		conf.println("</ideogram>");
		conf.close();
				
		conf=archiveFactory.openWriter(createFilename("ticks.conf"));
		conf.println("");
		conf.println("show_ticks        = no");
		conf.println("show_tick_labels  = no");
		conf.println("show_grid         = no");
		conf.println("");
		conf.println("<ticks>");
		conf.println("tick_label_font  = light");
		conf.println("radius           = dims(ideogram,radius_outer) + 180p");
		conf.println("label_offset     = 5p");
		conf.println("label_size       = 16p");
		conf.println("multiplier       = 1e-6");
		conf.println("color            = black");
		conf.println("thickness        = 1p");
		conf.println("");
		conf.println("# 25 Mb ticks, all chromosomes, with labels");
		conf.println("");
		conf.println("<tick>");
		conf.println("spacing        = 25u");
		conf.println("size           = 12p");
		conf.println("show_label     = yes");
		conf.println("format         = %d");
		conf.println("</tick>");
		conf.println("");
		conf.println("# 5 Mb ticks, all chromosomes, with labels");
		conf.println("# labels must be separated by at least 1px, which");
		conf.println("#   avoids overlap on human chrs");
		conf.println("<tick>");
		conf.println("label_separation = 1p");
		conf.println("spacing          = 5u");
		conf.println("size             = 7p");
		conf.println("show_label       = yes");
		conf.println("format           = %d");
		conf.println("</tick>");
		conf.println("");
		conf.println("# 5 Mb ticks with grid on rn1 and mm1, drawn only");
		conf.println("# for their grid (note size=0p)");
		conf.println("<tick>");
		conf.println("chromosomes_display_default = no");
		conf.println("chromosomes    = rn1;mm1");
		conf.println("spacing        = 5u");
		conf.println("size           = 0p");
		conf.println("force_display  = yes");
		conf.println("grid_start     = 0.45r");
		conf.println("grid_end       = dims(ideogram,radius_outer) + 180p");
		conf.println("grid_color     = grey");
		conf.println("grid_thickness = 1p");
		conf.println("grid           = yes");
		conf.println("</tick>");
		conf.println("");
		conf.println("# 20% relative ticks on human chromosomes");
		conf.println("");
		conf.println("<tick>");
		conf.println("chromosomes    = -rn1;-mm1");
		conf.println("radius         = 0.95r");
		conf.println("spacing_type   = relative");
		conf.println("rspacing       = 0.20");
		conf.println("size           = 6p");
		conf.println("show_label     = yes");
		conf.println("label_relative = yes");
		conf.println("rmultiplier    = 100");
		conf.println("format         = %d");
		conf.println("suffix         = %");
		conf.println("");
		conf.println("skip_last_label= yes");
		conf.println("");
		conf.println("grid_start     = 0.885r");
		conf.println("grid_end       = 0.95r");
		conf.println("grid_color     = grey");
		conf.println("grid_thickness = 1p");
		conf.println("grid           = yes");
		conf.println("</tick>");
		conf.println("");
		conf.println("# relative ticks at start and end of human chromosomes");
		conf.println("# with grid");
		conf.println("");
		conf.println("<tick>");
		conf.println("spacing_type   = relative");
		conf.println("rspacing       = 1");
		conf.println("size           = 6p");
		conf.println("grid_start     = 0.45r");
		conf.println("grid_end       = dims(ideogram,radius_outer) + 180p");
		conf.println("grid_color     = grey");
		conf.println("grid_thickness = 1p");
		conf.println("grid           = yes");
		conf.println("</tick>");
		conf.println("");
		conf.println("<tick>");
		conf.println("chromosomes    = -rn1;-mm1");
		conf.println("radius         = 0.95r");
		conf.println("spacing_type   = relative");
		conf.println("rspacing       = 0.10");
		conf.println("size           = 3p");
		conf.println("show_label     = no");
		conf.println("");
		conf.println("grid_start     = 0.885r");
		conf.println("grid_end       = 0.95r");
		conf.println("grid_color     = lgrey");
		conf.println("grid_thickness = 1p");
		conf.println("grid           = yes");
		conf.println("</tick>");
		conf.println("");
		conf.println("# 25% relative ticks on human chromosomes");
		conf.println("");
		conf.println("<tick>");
		conf.println("chromosomes    = -rn1;-mm1");
		conf.println("radius         = 0.82r");
		conf.println("spacing_type   = relative");
		conf.println("rspacing       = 0.25");
		conf.println("size           = 6p");
		conf.println("show_label     = yes");
		conf.println("label_relative = yes");
		conf.println("rmultiplier    = 100");
		conf.println("format         = %d");
		conf.println("");
		conf.println("skip_last_label= yes");
		conf.println("");
		conf.println("grid_start     = 0.755r");
		conf.println("grid_end       = 0.82r");
		conf.println("grid_color     = grey");
		conf.println("grid_thickness = 1p");
		conf.println("grid           = yes");
		conf.println("</tick>");
		conf.println("");
		conf.println("</ticks>");
		conf.close();
		
		conf = archiveFactory.openWriter(createFilename("circos.conf"));
		conf.println("karyotype="+createFilename("karyotype.txt"));
		conf.println("chromosomes_units=1000000");
		conf.println("<<include "+createFilename("ideogram.conf")+">>");
		conf.println("<<include "+createFilename("ticks.conf")+">>");

		
		final Map<String,SampleWriters> sample2writers=new HashMap<>(header.getNGenotypeSamples());
		for(final String sampleName: header.getSampleNamesInOrder())
			{
			sample2writers.put(sampleName, new SampleWriters(sampleName));
			}
		
		while(r.hasNext())
			{
			final VariantContext vc = r.next();
			for(final SampleWriters sw:sample2writers.values())
				{
				sw.visit(vc);
				}
			}
		r.close();r=null;
		LOG.info("writing conf");
		final Set<String> seen_chromosomes =new HashSet<>(dict.size());
		for(final SampleWriters sw:sample2writers.values())
			{
			seen_chromosomes.addAll(sw.bdns.stream().map(SV->SV.contig).collect(Collectors.toSet()));
			seen_chromosomes.addAll(sw.bdns.stream().map(SV->SV.bnb_contig).collect(Collectors.toSet()));
			seen_chromosomes.addAll(sw.duplications.stream().map(SV->SV.contig).collect(Collectors.toSet()));
			seen_chromosomes.addAll(sw.invertions.stream().map(SV->SV.contig).collect(Collectors.toSet()));
			seen_chromosomes.addAll(sw.deletions.stream().map(SV->SV.contig).collect(Collectors.toSet()));
			}
		
		pw = archiveFactory.openWriter(createFilename("karyotype.txt"));
		for(int i=0;i< dict.size();++i)
			{
			final SAMSequenceRecord ssr = dict.getSequence(i);
			if(!seen_chromosomes.contains(ssr.getSequenceName())) continue;
			pw.print("chr - ");
			pw.print(ssr.getSequenceName());
			pw.print(" ");
			pw.print(ssr.getSequenceName());
			pw.print(" 0 ");
			pw.print(ssr.getSequenceLength());
			pw.print(" ");
			pw.println(i%2==0?"blue":"grey");
			}
		pw.flush();pw.close();

		
		
		Function<String,String> sample2color = S->{
			int i= header.getSampleNamesInOrder().indexOf(S);
			if(i!=-1)
				{
				if(i%2==0) i=(( header.getSampleNamesInOrder().size()-1)-i);
				
				int g = 30+ (int)((i/(double)header.getSampleNamesInOrder().size())*50.0);
				
				return "gray"+(int)g;
				}
			
			
			
			return (i%2==0?"green":"red");
			};
		final double minR=0.1;
		final double maxR=0.95;
		final double numberOfTracks= sample2writers.values().stream().mapToInt(X->X.getNumberOfTracks()).sum();
		//final double maxSu = sample2writers.values().stream().mapToInt(X->X.getMaxSupportingEvidences()).max().orElse(1);
		double radius=minR;
		final double dr=(maxR-minR)/numberOfTracks;
		
		LOG.info("N tracks = "+numberOfTracks+ " dr="+dr);

		conf.println("<plots>");
		int maxSu =  sample2writers.values().stream().
				flatMap(S->S.deletions.stream()).
				mapToInt(SV->SV.su).max().orElse(0);
		for(final SampleWriters sw:sample2writers.values())
			{
			if(maxSu==0 || sw.deletions.isEmpty()) continue;
			String file =  sw.getDeletionFilename();
			pw = archiveFactory.openWriter(file);
			for(final SVDel sv: sw.deletions) sv.print(pw);
			pw.flush();
			pw.close();
			
			conf.println("<plot>");
			conf.println("type = histogram");
			conf.println("file = " + file);
			conf.println("min = 0");
			conf.println("max = " + maxSu);
			conf.println("orientation = out");
			conf.println("layers = 1");
			conf.println("margin = 0.02u");
			conf.println("color = black_a4");
			conf.println("r0 = "+radius+"r");
			conf.println("r1 = "+(radius+dr)+"r");
			conf.println("fill_color = " +sample2color.apply(sw.sampleName));
			conf.println("stroke_color = grey");// +sample2color.apply(sw.sampleName));
			conf.println("thickness = 1");
			conf.println("extend_bin  = no");
			//conf.println("<backgrounds>\n<background>\ncolor ="+(zz%2==0?"pink":"yellow")+"\n</background>\n</backgrounds>");
			//conf.println("<axes>\n<axis>\nspacing   = "+(dr/10.0)+"\ncolor=lgrey\nthickness=2\n</axis>\n</axes>");


			
			conf.println("</plot>");
			radius+=dr;
			}
		
		maxSu =  sample2writers.values().stream().
				flatMap(S->S.invertions.stream()).
				mapToInt(SV->SV.su).max().orElse(0);
		
		for(final SampleWriters sw:sample2writers.values())
			{
			if(maxSu==0 || sw.invertions.isEmpty()) continue;
			String file =  sw.getInvertionFilename();
			pw = archiveFactory.openWriter(file);
			for(final SVInv sv: sw.invertions) sv.print(pw);
			pw.flush();
			pw.close();
			
			conf.println("<plot>");
			conf.println("type = histogram");
			conf.println("file = " + file);
			conf.println("min = 0");
			conf.println("max = " + maxSu);
			conf.println("color = black_a4");
			conf.println("r0 = "+radius+"r");
			conf.println("r1 = "+(radius+dr)+"r");
			conf.println("fill_color = " +sample2color.apply(sw.sampleName));
			conf.println("stroke_color = " +sample2color.apply(sw.sampleName));
			conf.println("thickness = 1");
			conf.println("extend_bin  = no");
			conf.println("</plot>");
			radius+=dr;
			}
		
		maxSu =  sample2writers.values().stream().
				flatMap(S->S.duplications.stream()).
				mapToInt(SV->SV.su).max().orElse(0);
		
		for(final SampleWriters sw:sample2writers.values())
			{
			if(maxSu==0 || sw.duplications.isEmpty()) continue;
			String file =  sw.getDuplicationFilename();
			pw = archiveFactory.openWriter(file);
			for(final SVDup sv: sw.duplications) sv.print(pw);
			pw.flush();
			pw.close();
			
			conf.println("<plot>");
			conf.println("type = histogram");
			conf.println("file = " + file);
			conf.println("min = 0");
			conf.println("max = " + maxSu);
			conf.println("color = black_a4");
			conf.println("r0 = "+radius+"r");
			conf.println("r1 = "+(radius+dr)+"r");
			conf.println("fill_color = " +sample2color.apply(sw.sampleName));
			conf.println("stroke_color = " +sample2color.apply(sw.sampleName));
			conf.println("thickness = 1");
			conf.println("extend_bin  = no");
			conf.println("</plot>");
			radius+=dr;
			}
		
		conf.println("</plots>");
		
		conf.println("<links>");
		for(final SampleWriters sw:sample2writers.values())
			{
			if(sw.bdns.isEmpty()) continue;
				
			pw = archiveFactory.openWriter(sw.getBnbFilename());
			for(final SVBnB sv: sw.bdns) sv.print(pw);
			pw.flush();
			pw.close();
			
			conf.println("<link>");
			conf.println("file="+sw.getBnbFilename());
			conf.println("radius="+radius+"r");
			conf.println("bezier_radius = "+(radius+dr)+"r");
			conf.println("fill_color = " +sample2color.apply(sw.sampleName));
			conf.println("stroke_color = " +sample2color.apply(sw.sampleName));
			conf.println("</link>");
			radius+=dr;
				
			}
		conf.println("</links>");
		
		
		sample2writers.clear();
		
		conf.println("");
		conf.println("<image>");
		conf.println("# Included from Circos distribution.");
		conf.println("<<include etc/image.conf>>");
		conf.println("</image>");
		conf.println("");
		conf.println("<<include etc/colors_fonts_patterns.conf>>");
		conf.println("");
		conf.println("# Debugging, I/O an dother system parameters");
		conf.println("# Included from Circos distribution.");
		conf.println("<<include etc/housekeeping.conf>>");
		conf.println("");
		conf.flush();conf.close();conf=null;
		
		archiveFactory.close();
		archiveFactory=null;
		LOG.info("done");
		return 0;
	} catch(final Exception err)
	{
		LOG.error(err);
		return -1;
	}
		finally {
			CloserUtil.close(conf);
			CloserUtil.close(archiveFactory);
			CloserUtil.close(r);
		
		}
	}
	
	
public static void main(String[] args) {
	new LumpyVcfToCircos().instanceMainWithExit(args);
}
}
