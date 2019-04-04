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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;


/**
BEGIN_DOC

## Example

```
$ cat groups.txt

G1	10:112583204-112583210
G2	11
G3	12:1234-1235 13:20-30
```


```
$ java -jar dist/splitvcf.jar  -o tmp__GROUPID__.vcf.gz -g groups.txt in.vcf
$ ls tmp*
tmpG1.vcf.gz
tmpG2.vcf.gz
tmpG3.vcf.gz
tmpOTHER.vcf.gz
```


END_DOC
 
 */
@Program(name="splitvcf",
description="split a vcf using a named list of intervals...",
keywords={"vcf"})
public class SplitVcf
	extends Launcher
	{
	private static final Logger LOG = Logger.build(SplitVcf.class).make();

	private final static String REPLACE_GROUPID="__GROUPID__";
	private final java.util.Map<String,SplitGroup> name2group=new java.util.HashMap<>();
	private final IntervalTreeMap<Set<SplitGroup>> interval2group = new IntervalTreeMap<>();
	private SplitGroup underminedGroup=null;

	
	
	@Parameter(names={"-g","--groupfile"},description=
			"Chromosome group file. Intervals are 1-based. "
			+ "If undefined, splitvcf will use the sequence dictionary to output one vcf per contig."
			)
	private File chromGroupFile = null;
	
	@Parameter(names={"-u","--unmapped"},description="unmapped interval name")
	private String UNDERTERMINED_NAME = "OTHER";
	
	@Parameter(names={"-m","--multi"},description="if set, allow one variant to be mapped on multiple chromosome group (the record is duplicated)")
	private boolean multiIntervalEnabled = false;

	@Parameter(names={"-o","--out"},description="Output filename. Name must contain '__GROUPID__' ",
			required=true)
	private File outputFile = null;
	
	
	/** spit group */
	private class SplitGroup implements Closeable
		{
		final String groupName;
		VCFHeader header=null;
		VariantContextWriter _writer;
		
		SplitGroup(final String groupName)
			{
			this.groupName=groupName;
			}
		

		@Override
		public void close() {
			if(_writer!=null) CloserUtil.close(_writer);
			_writer=null;
			
		}
		@Override
		public int hashCode() {
			return this.groupName.hashCode();
			}
		
		@Override
		public boolean equals(Object obj) {
			if(obj==null) return false;
			if(obj == this) return true;
			return  this.groupName.equals(((SplitGroup)obj).groupName);
			}

		
		public File getFile()
			{
			return new File(
					SplitVcf.this.outputFile.getParentFile(),
					SplitVcf.this.outputFile.getName().replaceAll(
							SplitVcf.REPLACE_GROUPID,
							this.groupName
							));
			}
		
		public void open(final VCFHeader src)
			{	
	
			final File fileout=getFile();
			LOG.info("opening VCF file \""+fileout+"\" for writing");
			final File parent=fileout.getParentFile();
			if(parent!=null) {
				parent.mkdirs();
			}
	
			this.header= new VCFHeader(src);
			this.header.addMetaDataLine(new VCFHeaderLine("SplitVcf.GroupName", this.groupName));
			try {
				this._writer = VCFUtils.createVariantContextWriter(fileout);
			} catch (IOException e) {
				throw new RuntimeIOException(e);
			}
			this._writer.writeHeader(this.header);
			}
		
		@Override
		public String toString() {
			return groupName;
			}
		}
	
	
		
	public SplitVcf()
		{
		}
	
	
	
	private Set<SplitGroup> getGroupsFromInterval(final Interval interval) {
		final Set<SplitGroup> groups= new HashSet<>();
		for(final Set<SplitGroup> G:this.interval2group.getOverlapping(interval))
			{
			groups.addAll(G);
			}
		if(groups.isEmpty()) return Collections.emptySet();
		if(!this.multiIntervalEnabled && groups.size()!=1) {
			throw new IllegalStateException("got two group for uniq interval "+interval+" ? "+groups);
		}
		return groups;
		}
	
	
	private void putInterval(final String groupName,final Interval interval)
		{
		Collection<SplitGroup> splitgroups = this.getGroupsFromInterval(interval);
		if(!this.multiIntervalEnabled)
			{
			for(final SplitGroup splitgroup:splitgroups)
				{
				if(!splitgroup.groupName.equals(groupName))
					{
					throw new IllegalArgumentException(
							"chrom "+interval+" already used in "+splitgroup.groupName+
							". use option to enable multiple groups/interval");
					}
				}
			}
		
		SplitGroup splitgroup = this.name2group.get(groupName);
			
		if(splitgroup==null)
			{
			splitgroup = new SplitGroup(groupName);
			this.name2group.put(groupName,splitgroup);
			}
		Set<SplitGroup> L = this.interval2group.get(interval);
		if(L==null) {
			L=new HashSet<>();
			L.add(splitgroup);
		}
		this.interval2group.put(interval, L);
		}
				
	private int doWork(final String inputName) {
		if (this.outputFile==null || !this.outputFile.getName().contains(REPLACE_GROUPID)) {
			throw new JvarkitException.UserError("Output file pattern undefined or doesn't contain " + REPLACE_GROUPID + " : "
					+ this.outputFile);
		}
		if (!(this.outputFile.getName().endsWith(".vcf") || this.outputFile.getName().endsWith(".vcf.gz"))) {
			throw new JvarkitException.UserError("output file must end with '.vcf' or '.vcf.gz'");
		}
		BufferedReader r=null;
		VCFIterator in =null;
		try 
			{
			in = openVCFIterator(inputName);
			final SAMSequenceDictionary samSequenceDictionary = in.getHeader().getSequenceDictionary();
			if(samSequenceDictionary==null) {
				throw new JvarkitException.VcfDictionaryMissing(inputName==null?"<input>":inputName);
			}
			
			
			this.underminedGroup = new SplitGroup(UNDERTERMINED_NAME);
			this.name2group.put(UNDERTERMINED_NAME, this.underminedGroup);
	
			if (this.chromGroupFile != null)
				{
				r = IOUtils.openFileForBufferedReading(this.chromGroupFile);
				String line;
				while ((line = r.readLine()) != null) {
					if (line.isEmpty() || line.startsWith("#"))
						continue;
					final String tokens[] = line.split("[ \t,]+");
					final String groupName = tokens[0].trim();
					if (groupName.isEmpty())
						throw new JvarkitException.UserError("Empty group name in " + line);
					if (this.UNDERTERMINED_NAME.equals(groupName))
						throw new JvarkitException.UserError("Group cannot be named " + UNDERTERMINED_NAME);
					if (this.name2group.containsKey(groupName))
						throw new JvarkitException.UserError("Group defined twice " + groupName);
					for (int i = 1; i < tokens.length; i++) {
						String sequence;
						int start;
						int end;
						final String segment = tokens[i].trim();
	
						if (segment.isEmpty())
							continue;
	
						int colon = segment.indexOf(':');
						if (colon == -1) {
							final SAMSequenceRecord ssr = samSequenceDictionary.getSequence(segment);
							if (ssr == null) {
								throw new JvarkitException.ContigNotFoundInDictionary(segment,samSequenceDictionary);
							}
							sequence = segment;
							start = 1;
							end = ssr.getSequenceLength();
						} else {
							int hyphen = segment.indexOf('-', colon);
							if (hyphen == -1)
								throw new JvarkitException.UserError("Bad segment:" + segment);
							sequence = segment.substring(0, colon);
							if (samSequenceDictionary.getSequence(sequence) == null)
								throw new JvarkitException.UserError("Unknown chromosome , not in dict " + segment);
	
							start = Integer.parseInt(segment.substring(colon + 1, hyphen));
							end = Integer.parseInt(segment.substring(hyphen + 1));
						}
	
						final Interval interval = new Interval(sequence, start, end);
						this.putInterval(groupName, interval);
					}
				}
				r.close();r=null;
			} else {
				LOG.info("creating default split interval");
				for (final SAMSequenceRecord seq : samSequenceDictionary.getSequences()) {
					final String groupName = seq.getSequenceName();
					final Interval interval = new Interval(groupName, 1, seq.getSequenceLength());
					this.putInterval(groupName, interval);
				}
			}
			
			/* open all output vcf */
			for (final SplitGroup g : this.name2group.values()) {
				g.open(in.getHeader());
			}

			/* loop over vcf variations */
			final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(samSequenceDictionary);
			while (in.hasNext()) {
				final VariantContext record = progress.watch(in.next());
				final Interval interval = new Interval(record.getContig(), record.getStart(), record.getEnd());
				Set<SplitGroup> splitGroups = this.getGroupsFromInterval(interval);
				
				
				if (splitGroups.isEmpty()){
					splitGroups = Collections.singleton(this.underminedGroup);
				}
			
				
				for(final SplitGroup splitGroup:splitGroups)
					{
					splitGroup._writer.add(record);
					}
			}
	
			progress.finish();
			in.close();
	
			for (final SplitGroup g : this.name2group.values()) {
				g.close();
			}
			return 0;
			}
		catch(final Exception err) {
			for (final SplitGroup g : this.name2group.values()) {
				CloserUtil.close(g);
				if(in!=null) g.getFile().delete();
			}
			return -1;
		} finally {
			for (final SplitGroup g : this.name2group.values()) {
				CloserUtil.close(g);
			}
			CloserUtil.close(r);
			CloserUtil.close(in);
			this.name2group.clear();
			this.interval2group.clear();
		}
	}
	
	
	@Override
	public int doWork(final List<String> args) {
		return doWork(oneFileOrNull(args));
		}
	
	public static void main(final String[] args)
		{
		new SplitVcf().instanceMainWithExit(args);
		}

	}
