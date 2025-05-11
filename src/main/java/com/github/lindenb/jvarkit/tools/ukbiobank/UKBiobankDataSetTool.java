package com.github.lindenb.jvarkit.tools.ukbiobank;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.iterator.EqualIterator;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SortingCollection;

public class UKBiobankDataSetTool extends Launcher {
	private static final Logger LOG = Logger.of(UKBiobankDataSetTool.class);
	private UKBiobankDataSet dataset;
	@Parameter(names={"-d","--dataset"},description="Dataset basename", required = true)
	private String base;
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	
	
	private final List<UKBiobankDataSet.DataDict> predicates = new ArrayList<>();
	private final Map<String,Integer> column2index = new HashMap<>();
	
	private static class Triple {
		String subject;
		int predicate_idx;
		String value;
		int compare1(Triple o) {
			return subject.compareTo(o.subject);
			}
		int compare2(Triple o) {
			int i = compare1(o);
			if(i!=0) return i;
			i = Integer.compare(predicate_idx, o.predicate_idx);
			if(i!=0) return i;
			return value.compareTo(o.value);
			}
		}
	private class TripleCodec extends AbstractDataCodec<Triple> {
		@Override
		public void encode(DataOutputStream dos, Triple t) throws IOException {
			dos.writeUTF(t.subject);
			dos.writeInt(t.predicate_idx);
			dos.writeUTF(t.value);
			}
		@Override
		public Triple decode(DataInputStream dis) throws IOException {
			Triple t = new Triple();
			try {
				t.subject = dis.readUTF();
				}
			catch(EOFException err) {
				return null;
				}
			t.predicate_idx = dis.readInt();
			t.value = dis.readUTF();
			return t;
			}
		@Override
		public TripleCodec clone() {
			return new TripleCodec();
			}
		}
	
	private void dumpChildren(List<String> args) {
		final Set<UKBiobankDataSet.Coding> set = new HashSet<>();
		for(final String str: args) {
			int dot = str.indexOf('.');
			if(dot==-1) throw new IllegalArgumentException();
			final String code= str.substring(dot+1).trim();
			UKBiobankDataSet.Coding root= this.dataset.findCodingsByName(str.substring(0, dot)).
				stream().
				filter(C->C.getCode().equals(code)).
				findFirst().
				orElse(null);
			if(root==null) throw new IllegalArgumentException("Cannot find term "+code);
			set.addAll(root.getAllDescendants());
 			}
		PrintStream out = System.out;
		for(UKBiobankDataSet.Coding c:set) {
			out.print(c.getCodingName());
			out.print("\t");
			out.print(c.getCode());
			out.print("\t");
			out.print(c.getMeaning());
			out.print("\t");
			out.print(c.getParentCode());
			out.println();
			}
		out.flush();
		}
	
	private void make(List<String> args) throws IOException {
		SortingCollection<Triple> sorting = null;
		try {
			sorting = SortingCollection.newInstance(
					Triple.class,
					new TripleCodec(),
					(A,B)->A.compare2(B),
					writingSortingCollection.getMaxRecordsInRam(),
					writingSortingCollection.getTmpPaths()
					);
			sorting.setDestructiveIteration(true);
			CharSplitter tab = CharSplitter.TAB;
			for(Path p:IOUtils.unrollPaths(args)) {
				try(BufferedReader br=IOUtils.openPathForBufferedReading(p)) {
					String line = br.readLine();
					String[] tokens= tab.split(line);
					if(tokens[0].equals("participant.eid")) throw new IOException("first column should be participant.eid not "+tokens[0]);
					final List<UKBiobankDataSet.DataDict> cols = new ArrayList<>(tokens.length);
					cols.add(null);
					for(int i=1;i< tokens.length;i++) {
						UKBiobankDataSet.DataDict dd = dataset.getDataDictByColumn(tokens[i]);
						if(!this.column2index.containsKey(dd.getColumnName())) {
							 this.column2index.put(dd.getColumnName(),this.predicates.size());
							 this.predicates.add(dd);
							}
						cols.add(dd);
						}
					while((line=br.readLine())!=null) {
						tokens= tab.split(line);
						String participant = tokens[0];
						for(int i=1;i< tokens.length;i++) {
							if(tokens[i].isEmpty()) continue;
							final Triple t = new Triple();
							t.subject  = participant;
							t.predicate_idx = this.column2index.get(cols.get(i).getColumnName());
							t.value = tokens[i];
							sorting.add(t);
							}
						}
					}
				}
			sorting.doneAdding();
			
			try(ArchiveFactory archive = ArchiveFactory.open(Paths.get("xxx"))) {
				PrintWriter tableW = archive.openWriter("table.tsv");
				tableW.print("participant.eid");
				for(int x=0;x<this.predicates.size();++x) {
					tableW.print("\t");
					tableW.print(this.predicates.get(x).getTitle());
					}
				tableW.println();
				
				try(CloseableIterator<Triple> it0 = sorting.iterator()) {
					EqualIterator<Triple> it = new EqualIterator<>(it0,(A,B)->A.compare1(B));
					while(it.hasNext()) {
						final List<Triple> array = it.next();
						final Triple first = array.get(0);
						
						tableW.print(first.subject);
						for(int x=0;x<this.predicates.size();++x) {
							final int final_x=x;
							
							final String value= array.stream().
									filter(T->T.predicate_idx==final_x).
									map(T->T.value).
									findFirst().
									orElse("");
							tableW.print("\t");
							tableW.print(this.predicates.get(x).format1(value));
							}
						tableW.println();
						}
					}
				tableW.flush();
				tableW.close();
				}
			}
		catch(IOException err) {
			throw err;
			}
		}
	
	@Override
	public int doWork(List<String> args) {
		try {
			this.dataset = UKBiobankDataSet.load(this.base);
			return 0;
			}
		catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		}

	public static void main(final String[] args) {
		new UKBiobankDataSetTool().instanceMainWithExit(args);

	}
}
