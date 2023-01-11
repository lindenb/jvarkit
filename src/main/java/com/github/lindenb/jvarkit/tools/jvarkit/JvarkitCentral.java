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
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.tools.bam2graphics.WGSCoveragePlotter;
import com.github.lindenb.jvarkit.tools.bamstats04.BamStats05;
import com.github.lindenb.jvarkit.tools.bioalcidae.BioAlcidaeJdk;
import com.github.lindenb.jvarkit.tools.biostar.Biostar154220;
import com.github.lindenb.jvarkit.tools.gnomad.VcfGnomad;
import com.github.lindenb.jvarkit.tools.misc.VcfHead;
import com.github.lindenb.jvarkit.tools.misc.VcfTail;
import com.github.lindenb.jvarkit.tools.vcf2table.VcfToTable;
import com.github.lindenb.jvarkit.tools.vcffilterjs.VcfFilterJdk;
import com.github.lindenb.jvarkit.tools.vcffilterso.VcfFilterSequenceOntology;
import com.github.lindenb.jvarkit.tools.vcfpolyx.VCFPolyX;
import com.github.lindenb.jvarkit.tools.vcfviewgui.SwingBamCov;
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
			return p==null?getName():p.description();
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
		}
	private final List<Command> commands = new ArrayList<>();
	private Command command(Class<?> clazz) {
		Command c = new Command(clazz);
		this.commands.add(c);
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
		final int maxlenname = this.commands.stream().
				filter(P->!P.isHidden()).
				mapToInt(S->S.getName().length()).max().orElse(1);
		final int maxlendesc = this.commands.stream().
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
		for(Command c: this.commands) {
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
		command(BamStats05.class);
		command(BioAlcidaeJdk.class);
		command(Biostar154220.class);
		command(VCFPolyX.class);
		command(VcfToTable.class);
		command(VcfGnomad.class);
		command(VcfFilterSequenceOntology.class);
		command(WGSCoveragePlotter.class);
		command(SwingVcfView.class);
		command(SwingBamCov.class);
		command(VcfFilterJdk.class);
		command(SwingVcfJexlFilter.class);
		command(VcfHead.class);
		command(VcfTail.class);
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
		Collections.sort(this.commands,(A,B)->A.getName().compareTo(B.getName()));
		final String userCmd = args[0];
		final Command cmd = commands.stream().
				filter(C->C.getName().equals(args[0])).
				findAny().
				orElse(null);
		if(cmd!=null) {
			final String[] args2=new String[args.length-1];
			System.arraycopy(args, 1, args2, 0, args2.length);
			cmd.execute(args2);
			}
		else
			{
			LOG.error("jvarkit command \""+userCmd+"\" not found.\n"
					+ "Available commands are: " + this.commands.stream().map(S->S.getName()).collect(Collectors.joining(", ")));
			System.exit(-1);
			}
		}
	
    public static void main(String[] args) {
        new JvarkitCentral().run(args);
    }
}
