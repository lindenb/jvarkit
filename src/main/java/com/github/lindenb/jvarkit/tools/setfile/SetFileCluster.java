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
package com.github.lindenb.jvarkit.tools.setfile;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.setfile.SetFileRecord;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;

/**
BEGIN_DOC


TODO



END_DOC

**/
@Program(name="setfilecluster",
description="Cluster records of setfiles into files containing a sum to basepaires close to 'x' bp",
creationDate="20210125",
modificationDate="20240724",
keywords={"setfile","bed"},
jvarkit_amalgamion = true,
menu="BED Manipulation"
)
public class SetFileCluster extends AbstractSetFileTool {
	private static final Logger LOG = Logger.of(SetFileCluster.class);

	@Parameter(names={"-S","--size"},description="number of bases max per bin. (or specify --jobs). "+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.LongStringConverter.class,splitter=NoSplitter.class)
	private long long_length_per_bin=-1L;
	@Parameter(names={"-J","--jobs"},description="number of clusters. (or specify --size)")
	private int number_of_jobs=-1;

	private static class Cluster {
		final List<SetFileRecord> records = new ArrayList<>();
		long sum_length = 0L;
		void add(SetFileRecord rec) {
			this.records.add(rec);
			this.sum_length += rec.getLongSumOfLengthOnReference();
		}
		long getSumLength(final SetFileRecord malus) {
			return this.sum_length + (malus==null?0:malus.getLongSumOfLengthOnReference());
		}

	}
	
	
	@Override
	public int doWork(final List<String> args) {
		
		try {
			if(this.outputFile==null) {
				LOG.equals("output file must be defined.");
				return -1;
				}
			if(this.number_of_jobs<1 && this.long_length_per_bin<1L) {
				LOG.error("at least --jobs or --size must be specified.");
				return -1;
				}
			if(this.number_of_jobs>0 &&  this.long_length_per_bin>0) {
				LOG.error(" --jobs OR --size must be specified. Not both.");
				return -1;
				}
			
			final List<Cluster> clusters = new ArrayList<>();
		
			try(CloseableIterator<SetFileRecord> iter = openSetFileIterator(args)) {
				final List<SetFileRecord> records = iter.stream().
						filter(R->!R.isEmpty()).
						sorted((A,B)->Long.compare(B.getLongSumOfLengthOnReference(), A.getLongSumOfLengthOnReference())).
					collect(Collectors.toCollection(LinkedList::new));
				while(!records.isEmpty()) {
					final SetFileRecord first = records.remove(0);
					if(number_of_jobs>0) {
						if(clusters.size() < this.number_of_jobs) {
							final Cluster c = new Cluster();
							c.add(first);
							}
						else {
							int best_idx=-1;
							double best_length=-1;
							for(int y=0;y< clusters.size();++y) {
								final double total_length = clusters.get(y).getSumLength(first);
								if(best_idx==-1 ||total_length<best_length ) {
									best_idx=y;
									best_length = total_length;
									}
								}
							clusters.get(best_idx).add(first);
							}
						}
					else {
						int y=0;
						while(y<clusters.size()) {
								final Cluster cluster = clusters.get(y);
								if(cluster.getSumLength(first)<=this.long_length_per_bin) {
									cluster.add(first);
									break;
									}
								y++;
								}
						if(y==clusters.size()) {
							final Cluster cluster = new Cluster();
							cluster.add(first);
							clusters.add(cluster);
							}
						}
					}// end wile !records.isEmpty
				}// end open
			int clusterid = 0;
			try(final ArchiveFactory archive = ArchiveFactory.open(this.outputFile)) {
				for(final Cluster cluster : clusters) {
					Collections.sort(cluster.records,(A,B)->{
						final Locatable s1  = A.get(0);
						final Locatable s2  = B.get(0);
						final int i = getSorter().compare(s1, s2);
						if(i!=0) return i;
						return A.getName().compareTo(B.getName());
						});
					
					final String filename = String.format("cluster.%05d"+SetFileRecord.FILE_EXTENSION,clusterid);
					try(PrintWriter pw = archive.openWriter(filename)) {
						for(final SetFileRecord rec: cluster.records) {
							print(pw,rec);
							}
						pw.flush();
						}
					LOG.info(filename+" "+cluster.getSumLength(null)+"bp");
					++clusterid;
					}
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	public static void main(final String[] args) {
		new SetFileCluster().instanceMainWithExit(args);
	}

}
