/*
The MIT License (MIT)

Copyright (c) 2023 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.jvarkit;

import java.io.PrintStream;
import java.lang.reflect.Method;
import java.util.Map;
import java.util.TreeMap;

import com.github.lindenb.jvarkit.tools.allelebalance.VcfAlleleBalance;
import com.github.lindenb.jvarkit.tools.backlocate.BackLocate;
import com.github.lindenb.jvarkit.tools.bam2graphics.Bam2Raster;
import com.github.lindenb.jvarkit.tools.bam2graphics.CoveragePlotter;
import com.github.lindenb.jvarkit.tools.bam2graphics.LowResBam2Raster;
import com.github.lindenb.jvarkit.tools.bam2graphics.WGSCoveragePlotter;
import com.github.lindenb.jvarkit.tools.bam2svg.BamToSVG;
import com.github.lindenb.jvarkit.tools.bam2svg.WesCnvSvg;
import com.github.lindenb.jvarkit.tools.bam2xml.Bam2Xml;
import com.github.lindenb.jvarkit.tools.bamstats04.BamStats05;
import com.github.lindenb.jvarkit.tools.bedtools.BedCluster;
import com.github.lindenb.jvarkit.tools.bedtools.BedNonOverlappingSet;
import com.github.lindenb.jvarkit.tools.bioalcidae.BioAlcidaeJdk;
import com.github.lindenb.jvarkit.tools.biostar.Biostar103303;
import com.github.lindenb.jvarkit.tools.biostar.Biostar130456;
import com.github.lindenb.jvarkit.tools.biostar.Biostar139647;
import com.github.lindenb.jvarkit.tools.biostar.Biostar145820;
import com.github.lindenb.jvarkit.tools.biostar.Biostar154220;
import com.github.lindenb.jvarkit.tools.biostar.Biostar165777;
import com.github.lindenb.jvarkit.tools.biostar.Biostar170742;
import com.github.lindenb.jvarkit.tools.biostar.Biostar172515;
import com.github.lindenb.jvarkit.tools.biostar.Biostar173114;
import com.github.lindenb.jvarkit.tools.biostar.Biostar175929;
import com.github.lindenb.jvarkit.tools.biostar.Biostar178713;
import com.github.lindenb.jvarkit.tools.biostar.Biostar214299;
import com.github.lindenb.jvarkit.tools.biostar.Biostar234081;
import com.github.lindenb.jvarkit.tools.biostar.Biostar234230;
import com.github.lindenb.jvarkit.tools.biostar.Biostar251649;
import com.github.lindenb.jvarkit.tools.biostar.Biostar322664;
import com.github.lindenb.jvarkit.tools.biostar.Biostar332826;
import com.github.lindenb.jvarkit.tools.biostar.Biostar336589;
import com.github.lindenb.jvarkit.tools.biostar.Biostar352930;
import com.github.lindenb.jvarkit.tools.biostar.Biostar398854;
import com.github.lindenb.jvarkit.tools.biostar.Biostar404363;
import com.github.lindenb.jvarkit.tools.biostar.Biostar480685;
import com.github.lindenb.jvarkit.tools.biostar.Biostar489074;
import com.github.lindenb.jvarkit.tools.biostar.Biostar497922;
import com.github.lindenb.jvarkit.tools.biostar.Biostar59647;
import com.github.lindenb.jvarkit.tools.biostar.Biostar76892;
import com.github.lindenb.jvarkit.tools.biostar.Biostar77288;
import com.github.lindenb.jvarkit.tools.biostar.Biostar77828;
import com.github.lindenb.jvarkit.tools.biostar.Biostar78285;
import com.github.lindenb.jvarkit.tools.biostar.Biostar81455;
import com.github.lindenb.jvarkit.tools.biostar.Biostar84452;
import com.github.lindenb.jvarkit.tools.biostar.Biostar84786;
import com.github.lindenb.jvarkit.tools.biostar.Biostar86363;
import com.github.lindenb.jvarkit.tools.biostar.Biostar86480;
import com.github.lindenb.jvarkit.tools.biostar.Biostar90204;
import com.github.lindenb.jvarkit.tools.biostar.Biostar9462889;
import com.github.lindenb.jvarkit.tools.biostar.Biostar9469733;
import com.github.lindenb.jvarkit.tools.biostar.Biostar9501110;
import com.github.lindenb.jvarkit.tools.dbsnp.BuildDbsnp;
import com.github.lindenb.jvarkit.tools.gnomad.VcfGnomad;
import com.github.lindenb.jvarkit.tools.gtf.GtfToBed;
import com.github.lindenb.jvarkit.tools.minibam.MakeMiniBam;
import com.github.lindenb.jvarkit.tools.misc.ConvertBamChromosomes;
import com.github.lindenb.jvarkit.tools.misc.ConvertBedChromosomes;
import com.github.lindenb.jvarkit.tools.misc.ConvertVcfChromosomes;
import com.github.lindenb.jvarkit.tools.misc.VcfHead;
import com.github.lindenb.jvarkit.tools.misc.VcfTail;
import com.github.lindenb.jvarkit.tools.pubmed.Pubmed404;
import com.github.lindenb.jvarkit.tools.pubmed.PubmedCodingLanguages;
import com.github.lindenb.jvarkit.tools.pubmed.PubmedGender;
import com.github.lindenb.jvarkit.tools.pubmed.PubmedGraph;
import com.github.lindenb.jvarkit.tools.sam2tsv.CnvTView;
import com.github.lindenb.jvarkit.tools.sam2tsv.PrettySam;
import com.github.lindenb.jvarkit.tools.sam2tsv.Sam2Tsv;
import com.github.lindenb.jvarkit.tools.vcf2table.VcfToTable;
import com.github.lindenb.jvarkit.tools.vcffilterjs.VcfFilterJdk;
import com.github.lindenb.jvarkit.tools.vcffilterso.VcfFilterSequenceOntology;
import com.github.lindenb.jvarkit.tools.vcfpolyx.VCFPolyX;
import com.github.lindenb.jvarkit.tools.vcfrebase.VcfRebase;
import com.github.lindenb.jvarkit.tools.vcfsplit.VcfSplitNVariants;
import com.github.lindenb.jvarkit.tools.vcfviewgui.SwingBamCov;
import com.github.lindenb.jvarkit.tools.vcfviewgui.SwingBamView;
import com.github.lindenb.jvarkit.tools.vcfviewgui.SwingVcfJexlFilter;
import com.github.lindenb.jvarkit.tools.vcfviewgui.SwingVcfView;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.StringUtil;

public class JvarkitCentral {
	private static final Logger LOG = Logger.build(JvarkitCentral.class).make();
	private class Command {
		private final Class<?> clazz;
		boolean hidden = false;
		Command(final Class<?> claz) {
			this.clazz = claz;
			}
		
		public boolean isHidden() {
			return hidden;
			}
		
		public String getName() {
			final Program p = this.clazz.getAnnotation(Program.class);
			return p==null?this.clazz.getSimpleName().toLowerCase():p.name();
			}
		public String getDescription() {
			final Program p = this.clazz.getAnnotation(Program.class);
			return p==null?getName():p.description().trim();
			}
		void execute(final String[] args) {
			try {
				final Method mainMethod = this.clazz.getMethod("main",String[].class);
				mainMethod.invoke(null, (Object)args);
				}
			catch(final Throwable err) {
				LOG.error(err);
				System.exit(-1);
				}
			}
		Command setHidden() {
			this.hidden = true;
			return this;
			}
		}
	private final Map<String,Command> commands = new TreeMap<>();
	private Command command(Class<?> clazz) {
		Command c = new Command(clazz);
		if(this.commands.containsKey(c.getName())) {
			return this.commands.get(c.getName());
			}
		this.commands.put(c.getName(),c);
		return c;
		}
	
	private void usage(PrintStream out) {
		out.println("JVARKIT");
		out.println("=======");
		out.println();
		out.println("Author      : Pierre Lindenbaum Phd. Institut du Thorax. Nantes. France.");
		out.println("Version     : " + JVarkitVersion.getInstance().getGitHash());
		out.println("Compilation : " + JVarkitVersion.getInstance().getCompilationDate());
		out.println("Github      : https://github.com/lindenb/jvarkit");
		out.println("Issues      : https://github.com/lindenb/jvarkit/issues");
		out.println();
		out.println("## Usage");
		out.println();
		out.println(" java -jar jvarkit.jar <command name> (other arguments)");
		out.println();
		out.println("## Options");
		out.println();
		out.println(" + --help show this screen");
		out.println(" + --version print version");
		out.println();
		out.println("## Tools");
		out.println();
		final int maxlenname = this.commands.values().stream().
				filter(P->!P.isHidden()).
				mapToInt(S->S.getName().length()).max().orElse(1);
		final int maxlendesc = this.commands.values().stream().
				filter(P->!P.isHidden()).
				mapToInt(S->S.getDescription().length()).max().orElse(1);

		out.print("| ");
		out.print(String.format("%"+maxlenname+"s","Name"));
		out.print(" | ");
		out.print(String.format("%-"+maxlendesc+"s","Description"));
		out.println(" |");
		out.print("|-");
		out.print(StringUtil.repeatCharNTimes('-', maxlenname));
		out.print("-|-");
		out.print(StringUtil.repeatCharNTimes('-', maxlendesc));
		out.print("-|");
		out.println();
		for(Command c: this.commands.values()) {
			if(c.isHidden()) continue;
			out.print("| ");
			out.print(String.format("%"+maxlenname+"s",c.getName()));
			out.print(" | ");
			out.print(String.format("%-"+maxlendesc+"s",c.getDescription()));
			out.print(" |");
			out.println();
			}
		out.println();
	}
	
	private void run(String[] args) {
		command(BamToSVG.class);
		command(Bam2Xml.class);
		command(Bam2Raster.class);
		command(BamStats05.class);
		command(BackLocate.class);
		command(BedCluster.class);
		command(BedNonOverlappingSet.class);
		command(BioAlcidaeJdk.class);
		command(Biostar103303.class);
		command(Biostar154220.class);
		command(Biostar130456.class);
		command(Biostar139647.class);
		command(Biostar145820.class);
		command(Biostar165777.class);
		command(Biostar170742.class);
		command(Biostar172515.class);
		command(Biostar173114.class);
		command(Biostar175929.class);
		command(Biostar178713.class);
		command(Biostar214299.class);
		command(Biostar234081.class);
		command(Biostar234230.class);
		command(Biostar251649.class);
		command(Biostar322664.class);
		command(Biostar332826.class);
		command(Biostar336589.class);
		command(Biostar352930.class);
		command(Biostar398854.class);
		command(Biostar404363.class);
		command(Biostar480685.class);
		command(Biostar489074.class);
		command(Biostar497922.class);
		command(Biostar59647.class);
		command(Biostar76892.class);
		command(Biostar77288.class);
		command(Biostar77828.class);
		command(Biostar78285.class);
		command(Biostar81455.class);
		command(Biostar84452.class);
		command(Biostar84786.class);
		command(Biostar86363.class);
		command(Biostar86480.class);
		command(Biostar90204.class);
		command(Biostar9462889.class);
		command(Biostar9469733.class);
		command(Biostar9501110.class);
		command(BuildDbsnp.class);
		command(CnvTView.class);
		command(ConvertBamChromosomes.class);
		command(ConvertBedChromosomes.class);
		command(ConvertVcfChromosomes.class).setHidden();

		command(CoveragePlotter.class);
		command(GtfToBed.class);
		command(LowResBam2Raster.class);
		command(MakeMiniBam.class);
		command(Pubmed404.class);
		command(PubmedCodingLanguages.class);
		command(PubmedGender.class);
		command(PubmedGraph.class);
		command(PrettySam.class);
		command(VCFPolyX.class);
		command(VcfToTable.class);
		command(VcfGnomad.class);
		command(VcfFilterSequenceOntology.class);
		command(Sam2Tsv.class);
		command(SwingVcfView.class);
		command(SwingBamCov.class);
		command(SwingBamView.class);
		command(VcfFilterJdk.class);
		command(VcfAlleleBalance.class);
		command(SwingVcfJexlFilter.class);
		command(VcfRebase.class);
		command(VcfHead.class);
		command(VcfTail.class);
		command(VcfTail.class);
		command(VcfSplitNVariants.class);
		command(WesCnvSvg.class);
		command(WGSCoveragePlotter.class);
		
		if(args.length==0) {
			usage(System.err);
			return;
		}
		if(args[0].equals("--version") || args[0].equals("-v")) {
			System.out.println(JVarkitVersion.getInstance().getGitHash());
			return;
		}
		if(args[0].equals("--help") || args[0].equals("-h")) {
			usage(System.out);
			return;
			}
		if(args[0].startsWith("-")) {
			LOG.error("exepected a sub-jvarkit-program. but got "+args[0]+". type java -jar jvarkit.jar --help for more information.");
			System.exit(-1);
			}
		final String userCmd = args[0];
		final Command cmd = this.commands.get(userCmd);
		if(cmd!=null) {
			final String[] args2=new String[args.length-1];
			System.arraycopy(args, 1, args2, 0, args2.length);
			cmd.execute(args2);
			}
		else
			{
			LOG.error("jvarkit command \""+userCmd+"\" not found.\n"
					+ "Available commands are: " + String.join(", ",this.commands.keySet()));
			System.exit(-1);
			}
		}
	
    public static void main(String[] args) {
        new JvarkitCentral().run(args);
    }
}
