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

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.lang.reflect.Method;
import java.util.Arrays;
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
import com.github.lindenb.jvarkit.tools.bam2wig.Bam2Wig;
import com.github.lindenb.jvarkit.tools.bam2xml.Bam2Xml;
import com.github.lindenb.jvarkit.tools.bamstats04.BamStats05;
import com.github.lindenb.jvarkit.tools.basecoverage.BaseCoverage;
import com.github.lindenb.jvarkit.tools.bedtools.BedCluster;
import com.github.lindenb.jvarkit.tools.bedtools.BedNonOverlappingSet;
import com.github.lindenb.jvarkit.tools.bioalcidae.BioAlcidaeJdk;
import com.github.lindenb.jvarkit.tools.biostar.Biostar103303;
import com.github.lindenb.jvarkit.tools.biostar.Biostar105754;
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
import com.github.lindenb.jvarkit.tools.groupbygene.GroupByGene;
import com.github.lindenb.jvarkit.tools.gtf.GtfToBed;
import com.github.lindenb.jvarkit.tools.gvcf.FindGVCFsBlocks;
import com.github.lindenb.jvarkit.tools.minibam.MakeMiniBam;
import com.github.lindenb.jvarkit.tools.misc.AddLinearIndexToBed;
import com.github.lindenb.jvarkit.tools.misc.ConvertBamChromosomes;
import com.github.lindenb.jvarkit.tools.misc.ConvertBedChromosomes;
import com.github.lindenb.jvarkit.tools.misc.ConvertVcfChromosomes;
import com.github.lindenb.jvarkit.tools.misc.VCFShuffle;
import com.github.lindenb.jvarkit.tools.misc.VcfHead;
import com.github.lindenb.jvarkit.tools.misc.VcfTail;
import com.github.lindenb.jvarkit.tools.pubmed.Pubmed404;
import com.github.lindenb.jvarkit.tools.pubmed.PubmedCodingLanguages;
import com.github.lindenb.jvarkit.tools.pubmed.PubmedGender;
import com.github.lindenb.jvarkit.tools.pubmed.PubmedGraph;
import com.github.lindenb.jvarkit.tools.sam2tsv.CnvTView;
import com.github.lindenb.jvarkit.tools.sam2tsv.PrettySam;
import com.github.lindenb.jvarkit.tools.sam2tsv.Sam2Tsv;
import com.github.lindenb.jvarkit.tools.samgrep.SamGrep;
import com.github.lindenb.jvarkit.tools.samrmdupnames.SamRemoveDuplicatedNames;
import com.github.lindenb.jvarkit.tools.structvar.CoverageMatrix;
import com.github.lindenb.jvarkit.tools.structvar.VcfStrechToSvg;
import com.github.lindenb.jvarkit.tools.structvar.indexcov.SwingIndexCov;
import com.github.lindenb.jvarkit.tools.uniprot.UniprotToSvg;
import com.github.lindenb.jvarkit.tools.vcf2table.VcfToTable;
import com.github.lindenb.jvarkit.tools.vcfbigwig.VCFBigWig;
import com.github.lindenb.jvarkit.tools.vcfbigwig.VcfBigBed;
import com.github.lindenb.jvarkit.tools.vcffilterjs.VcfFilterJdk;
import com.github.lindenb.jvarkit.tools.vcffilterso.VcfFilterSequenceOntology;
import com.github.lindenb.jvarkit.tools.vcfpar.VcfPseudoAutosomalRegion;
import com.github.lindenb.jvarkit.tools.vcfpolyx.VCFPolyX;
import com.github.lindenb.jvarkit.tools.vcfrebase.VcfRebase;
import com.github.lindenb.jvarkit.tools.vcfsplit.VcfSplitNVariants;
import com.github.lindenb.jvarkit.tools.vcfsplitgene.VcfGeneSplitter;
import com.github.lindenb.jvarkit.tools.vcfviewgui.SwingBamCov;
import com.github.lindenb.jvarkit.tools.vcfviewgui.SwingBamView;
import com.github.lindenb.jvarkit.tools.vcfviewgui.SwingVcfJexlFilter;
import com.github.lindenb.jvarkit.tools.vcfviewgui.SwingVcfView;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StringUtil;

public class JvarkitCentral {
	private static final Logger LOG = Logger.build(JvarkitCentral.class).make();
	private class Command {
		private final Class<?> clazz;
		boolean hidden = false;
		Command(final Class<?> claz) {
			this.clazz = claz;
			}
		private Program program() {
				return this.clazz.getAnnotation(Program.class);
			}
		public boolean isHidden() {
			final Program p = program();
			if(p!=null && p.jvarkit_hidden()) return true;
			return hidden;
			}
		
		public String getName() {
			final Program p = program();
			return p==null?this.clazz.getSimpleName().toLowerCase():p.name();
			}
		public String getDescription() {
			final Program p = program();
			return p==null?getName():p.description().trim();
			}
		
		public String getCreationDate() {
			final Program p = program();
			return p==null?".":p.creationDate();
			}
		public String getModificationDate() {
			final Program p = program();
			return p==null?".":p.modificationDate();
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
		@SuppressWarnings("unused")
		Command setHidden() {
			this.hidden = true;
			return this;
			}
		
		
		
		void writeDoc(final File dir) throws IOException {
			final PrintStream old = System.out;
			final File filename = new File(dir,this.clazz.getSimpleName()+".md");
			final String[] cmd = new String[]{"--help","--helpFormat","markdown"};
			LOG.info((filename.exists()?"over":"")+"writing doc for "+getName()+" into "+filename);
			try (PrintStream out = new PrintStream(filename)) {
				System.setOut(out);
				if(Launcher.class.isAssignableFrom(this.clazz)) {
					Launcher instance= Launcher.class.cast(this.clazz.getConstructor().newInstance());
					final int  ret= instance.instanceMain(Arrays.asList(cmd));
					instance = null;
					if(ret!=0) {
						LOG.warn("return status!=0");
						}
					}
				else
					{
					execute(cmd);
					}
				}
			catch(Throwable err) {
				LOG.error("Cannot save doc for "+getName(), err);
				}
			finally
				{
				System.setOut(old);
				}
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
	
	
	private void usageMain(final PrintStream out) {
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
		out.println("  java -jar jvarkit.jar [options]");
		out.println("or");
		out.println("  java -jar jvarkit.jar <command name> (other arguments)");
		out.println();
		out.println("## Options");
		out.println();
		out.println(" + --help show this screen");
		out.println(" + --help-all show all commands, including the private ones.");
		out.println(" + --version print version");
		out.println();
		}
	
	private void usage(PrintStream out,boolean showHidden) {
		usageMain(out);
		
		out.println("## Tools");
		out.println();
		final int maxlenname = this.commands.values().stream().
				filter(P->showHidden || !P.isHidden()).
				mapToInt(S->S.getName().length()).max().orElse(1);
		final int maxlendesc = this.commands.values().stream().
				filter(P->showHidden || !P.isHidden()).
				mapToInt(S->S.getDescription().length()).max().orElse(1);

		out.print("| ");
		out.print(String.format("%"+maxlenname+"s","Name"));
		out.print(" | ");
		out.print(String.format("%-"+maxlendesc+"s","Description"));
		out.println(" |");
		out.print("|-");
		out.print(StringUtil.repeatCharNTimes('-', maxlenname));
		out.print("-+-");
		out.print(StringUtil.repeatCharNTimes('-', maxlendesc));
		out.print("-|");
		out.println();
		for(Command c: this.commands.values()) {
			if(!showHidden && c.isHidden()) continue;
			out.print("| ");
			out.print(String.format("%"+maxlenname+"s",c.getName()));
			out.print(" | ");
			out.print(String.format("%-"+maxlendesc+"s",c.getDescription()));
			out.print(" |");
			out.println();
			}
		out.println();
	}
	
	private String hyperlink(final String url) {
	return "["+url+"]("+url+")";
	}
	
	void writeReadTheDoc(final File dir) throws IOException {
		IOUtil.assertDirectoryIsWritable(dir);
		final File docs = new File(dir,"docs");
		IOUtil.assertDirectoryIsWritable(docs);
		try (PrintStream out = new PrintStream(new File(docs,"index.md"))) {
			usageMain(out);
			out.println("## Compilation Installation\n");
			out.println("Please, read [how to run and install jvarkit](JvarkitCentral.md)\n");
			out.println("## Tools\n");
			out.println("| Tool | Description | Creation | Update |");
			out.println("| ---: | :---------- | :------: | :----: |");
			for(Command c: this.commands.values()) {
				if(c.isHidden()) continue;
				out.print("| ");
				out.print(c.getName());
				out.print(" | ");
				out.print(c.getDescription());
				out.print(" | ");
				out.print(c.getCreationDate());
				out.print(" | ");
				out.print(c.getModificationDate());
				out.print(" |");
				out.println();
				}
			out.flush();
			}
		try (PrintStream out = new PrintStream(new File(docs,"JvarkitCentral.md"))) {
			usageMain(out);
			
			out.append("## Compilation\n");
			out.append("\n");
			out.append("### Requirements / Dependencies\n");
			out.append("\n");
			out.append("* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )\n");
			out.append("\n");
			out.append("\n");
			out.append("### Download and Compile\n");
			out.append("\n");
			out.append("```bash\n");
			out.append("$ git clone \"https://github.com/lindenb/jvarkit.git\"\n");
			out.append("$ cd jvarkit\n");
			out.append("$ ./gradlew jvarkit\n");
			out.append("```\n");
			out.append("\n");
			out.append("The java jar file will be installed in the `dist` directory.\n");
			out.append("\n");
			
			out.append("## Contribute\n");
			out.append("\n");
			out.append("- Issue Tracker: "+hyperlink("http://github.com/lindenb/jvarkit/issues")+"\n");
			out.append("- Source Code: "+hyperlink("http://github.com/lindenb/jvarkit")+"\n");
			out.append("\n");
			out.append("## License\n");
			out.append("\n");
			out.append("The project is licensed under the MIT license.\n");
			out.append("\n");
			out.append("## Citing\n");
			out.append("\n");
			out.append("Should you cite **jvarkit** ? "+hyperlink("https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md")+"\n");
			out.append("\n");
			out.append("The current reference is:\n");
			out.append("\n");
			out.append(hyperlink("http://dx.doi.org/10.6084/m9.figshare.1425030")+"\n");
			out.append("\n");
			out.append("> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.\n");
			out.append("> "+hyperlink("http://dx.doi.org/10.6084/m9.figshare.1425030")+"\n");
			out.append("\n");

			out.flush();
			}
		try (PrintStream out = new PrintStream(new File(dir,"mkdocs.yml"))) {
				out.println("repo_url: \"https://github.com/lindenb/jvarkit\"");
				out.println("repo_name: \"GitHub\"");
				out.println("docs_dir: docs");
				out.println("");
				out.println("theme:");
				out.println("    name: readthedocs");
				out.println("    highlightjs: true");
				out.println("    html_show_sourcelink: false");
				out.println("    prev_next_buttons_location: both");
				out.println("");
				out.println("nav:");
				for(Command c:this.commands.values()) {
					out.println("        - \""+ c.getName() +"\" : "+ c.clazz.getSimpleName() +".md");
					//out.println("        - "VcfTail": VcfTail.md");
					}
				out.println();
				out.flush();
			}
		for(Command c:this.commands.values()) {
			try {
				c.writeDoc(docs);
				}
			catch(IOException err) {
				LOG.warn(err);
				}
			}
		}
	
	private void run(String[] args) {
		command(AddLinearIndexToBed.class);
		command(BackLocate.class);
		command(BamToSVG.class);
		command(Bam2Xml.class);
		command(Bam2Wig.class);
		command(Bam2Raster.class);
		command(BamStats05.class);
		command(BaseCoverage.class);
		command(BedCluster.class);
		command(BedNonOverlappingSet.class);
		command(BioAlcidaeJdk.class);
		command(Biostar103303.class);
		command(Biostar105754.class);
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
		command(CoverageMatrix.class);
		command(BuildDbsnp.class);
		command(CnvTView.class);
		command(ConvertBamChromosomes.class);
		command(ConvertBedChromosomes.class);
		command(ConvertVcfChromosomes.class);
		command(CoveragePlotter.class);
		command(FindGVCFsBlocks.class);
		command(GtfToBed.class);
		command(GroupByGene.class);
		command(LowResBam2Raster.class);
		command(MakeMiniBam.class);
		command(Pubmed404.class);
		command(PubmedCodingLanguages.class);
		command(PubmedGender.class);
		command(PubmedGraph.class);
		command(PrettySam.class);
		command(VCFBigWig.class);
		command(VcfBigBed.class);
		command(VCFPolyX.class);
		command(VcfToTable.class);
		command(VcfGnomad.class);
		command(VcfFilterSequenceOntology.class);
		command(Sam2Tsv.class);
		command(SamGrep.class);
		command(SamRemoveDuplicatedNames.class);
		command(SwingVcfView.class);
		command(SwingBamCov.class);
		command(SwingBamView.class);
		command(SwingIndexCov.class);
		command(SwingVcfJexlFilter.class);
		command(UniprotToSvg.class);
		command(VcfFilterJdk.class);
		command(VcfAlleleBalance.class);
		command(VcfPseudoAutosomalRegion.class);
		command(VcfRebase.class);
		command(VcfHead.class);
		command(VcfTail.class);
		command(VcfGeneSplitter.class);
		command(VcfSplitNVariants.class);
		command(VcfStrechToSvg.class);
		command(VCFShuffle.class);
		command(WesCnvSvg.class);
		command(WGSCoveragePlotter.class);
		
		if(args.length==0) {
			usage(System.err,false);
			return;
		}
		else if(args[0].equals("--version") || args[0].equals("-v")) {
			System.out.println(JVarkitVersion.getInstance().getGitHash());
			return;
		}
		else if(args[0].equals("--help") || args[0].equals("-h")) {
			usage(System.out, false);
			return;
			}
		else if(args[0].equals("--help-all")) {
			usage(System.out, true);
			return;
			}
		else if(args.length==2 && args[0].equals("--generate-doc")) {
			final File dir = new File(args[1]);
			IOUtil.assertDirectoryIsWritable(dir);
			for(Command c:this.commands.values()) {
				try {
					c.writeDoc(dir);
					}
				catch(IOException err) {
					LOG.warn(err);
					}
				}
			return;
			}
		else if(args.length==2 && args[0].equals("--readthedocs")) {
			final File dir = new File(args[1]);
			try {
				writeReadTheDoc(dir);
				}
			catch(IOException err) {
				LOG.error(err);
				System.exit(-1);
				}
			return;
			}
		else if(args[0].startsWith("-")) {
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
