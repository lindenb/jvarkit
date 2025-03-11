package com.github.lindenb.jvarkit.tools.regenie;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;

public class RegenieMakeSetFile extends Launcher {
	private static final Logger LOG = Logger.build(RegenieMakeSetFile.class).make();

	@Parameter(names={"-o","--output"},description="output dir.",required = true)
	private Path outputDir;

private class Output implements Closeable{
	PrintWriter annot;
	PrintWriter setfile;
	PrintWriter mask;
	final Set<String> prediction = new HashSet<>();
	Output(int index) throws IOException {
		final String prefix = String.format("chunk%03d.",index);
		this.annot = IOUtils.openPathForPrintWriter(outputDir.resolve(prefix+"annot.txt"));
		this.setfile = IOUtils.openPathForPrintWriter(outputDir.resolve(prefix+"setfile.txt"));
		this.mask = IOUtils.openPathForPrintWriter(outputDir.resolve(prefix+"mask.txt"));
		}
	@Override
	public void close() {
		this.annot.flush();
		this.setfile.flush();
		this.mask.flush();
		this.annot.close();
		this.setfile.close();
		this.mask.close();
		}
	}
private static class Target {
	String contig;
	int start;
	int end;
	int output_index=-1;
	}



@Override
public int doWork(List<String> args) {
	try {
		if(args.size()!=3) {
			LOG.error("expected: annot setfile mask");
			return -1;
			}
		final Pattern spaces_regex = Pattern.compile("[ \t]");
		final Map<String, Target> target_hash = new HashMap<>(50_000);
		final Path annotfilePath = Paths.get(args.get(0));
		final Path setfilePath = Paths.get(args.get(1));
		final Path maskPath = Paths.get(args.get(2));
		IOUtil.assertFileIsReadable(annotfilePath);
		IOUtil.assertFileIsReadable(setfilePath);
		IOUtil.assertFileIsReadable(maskPath);
		IOUtil.assertDirectoryIsWritable(this.outputDir);
		try(BufferedReader br=IOUtils.openPathForBufferedReading(annotfilePath)) {
			String line;
			while((line=br.readLine())!=null) {
				final String[] tokens = spaces_regex.split(line);
				final String id = tokens[0];
				final String contig = CharSplitter.COLON.split(id)[0];
				final String gene = tokens[1];
				final String gene_id = contig+"~"+gene;
				final int pos = Integer.parseInt(tokens[0]);
				Target t = target_hash.get(gene_id);
				if(t==null) {
					t =  new Target();
					t.contig = contig;
					t.start = pos;
					t.end = pos;
					target_hash.put(gene_id, t);
					}
				t.start = Math.min(t.start, pos);
				t.end = Math.max(t.end, pos);
				}
			}
		int max_index=-1;
		final IntervalTreeMap<List<Target>> treeMap = new IntervalTreeMap<>();
		for(Target t:target_hash.values()) {
			final Interval loc = new Interval(t.contig, t.start, t.end);
			final Set<Integer> indexes= treeMap.getOverlapping(loc).stream().
					flatMap(L->L.stream()).
					map(T->T.output_index).
					collect(Collectors.toSet());
			t.output_index = 0;
			while(indexes.contains(t.output_index)) {
				t.output_index++;
				}
			max_index=Math.max(max_index, t.output_index);
			List<Target> L = treeMap.get(loc);
			if(L==null) {
				L = new ArrayList<>();
				treeMap.put(loc, L);
				}	
			L.add(t);
			}
		for(int i=0;i<= max_index;++i) {
			try(Output output=new Output(i)) {
				try(BufferedReader br=IOUtils.openPathForBufferedReading(annotfilePath)) {
					String line;
					while((line=br.readLine())!=null) {
						final String[] tokens = spaces_regex.split(line);
						final String id = tokens[0];
						final String gene = tokens[1];
						final String gene_id = CharSplitter.COLON.split(id)[0]+"~"+gene;
						final Target t = target_hash.get(gene_id);
						if(t==null) throw new IllegalStateException(line+"\n"+gene_id);
						if( t.output_index!=i) continue;
						output.annot.println(line);
						output.prediction.add(tokens[2]);
						}
					}
				try(BufferedReader br=IOUtils.openPathForBufferedReading(setfilePath)) {
					String line;
					while((line=br.readLine())!=null) {
						final String[] tokens = spaces_regex.split(line);
						final String gene_id = tokens[1] + "~" + tokens[0];
						final Target t = target_hash.get(gene_id);
						if(t==null) throw new IllegalStateException(line+"\n"+gene_id);
						if(t.output_index!=i) continue;
						output.setfile.println(line);
						}
					}
				try(BufferedReader br=IOUtils.openPathForBufferedReading(maskPath)) {
					String line;
					while((line=br.readLine())!=null) {
						final String[] tokens = spaces_regex.split(line);
						final Set<String> masks = new HashSet<>(CharSplitter.COMMA.splitAsStringList(tokens[1]));
						masks.retainAll(output.prediction);
						if(!masks.isEmpty()) {
							output.mask.println(line);
							}
						}
					}
				}
			}
		return 0;
		}
	catch(Throwable err) {
		err.printStackTrace();
		return -1;
		}
	}

public static void main(String[] args) {
	new RegenieMakeSetFile().instanceMainWithExit(args);
	}


}
