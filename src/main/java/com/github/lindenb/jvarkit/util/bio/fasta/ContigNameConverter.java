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


History:
* 2017 creation

*/
package com.github.lindenb.jvarkit.util.bio.fasta;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.cram.build.CramIO;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;

public abstract class ContigNameConverter implements  Function<String, String> {
public static final String OPT_ON_NT_FOUND_DESC="Contig converter. I will do my best to convert the contig names (e.g 'chr1' -> '1'): But what should I do when comparing two dictionaries with different notations";	
protected abstract String find(final String contig);

public String getName() {
	return "Untitled";
}



/** return converted contig, or null */
@Override
public String apply(final String contig) {
	final String dest=find(contig);
	return dest;
	}		

public SAMSequenceRecord convertSAMSequenceRecord(final SAMSequenceRecord ssr) {
	final String newName=apply(ssr.getSequenceName());
	if(newName==null)
		{
		return null;
		}
	final SAMSequenceRecord ssr2=new SAMSequenceRecord(newName, ssr.getSequenceLength());
		
	return ssr2;
	}


public SAMSequenceDictionary convertDictionary(final SAMSequenceDictionary dict) {
	return new SAMSequenceDictionary(dict.getSequences().
		stream().
		map(S->convertSAMSequenceRecord(S)).
		filter(S->S!=null).
		collect(Collectors.toList())
		);
	}

public static  ContigNameConverter getIdentity() {
	return new ContigNameConverter()
			{
			@Override
			protected String find(final String contig) {
				return contig;
				}
			@Override
			public String getName() {
				return "Identity";
			}
			};
	}

private static class MapBasedContigNameConverter extends ContigNameConverter
	{
	protected String name="Untitled";
	protected final Map<String,String> map;
	
	MapBasedContigNameConverter() {
		this.map=new HashMap<>();
		}
	void put(final String src,final String dest) {
		this.map.put(src, dest);
		this.map.put(dest, dest);
		}
	@Override
	protected String find(final String contig) {
		return this.map.getOrDefault(contig,null);
		}
	@Override
	public String getName() {
		return name;
		}
	
	}

public static final String OPT_DICT_OR_MAPPING_FILE_DESC="Chromosome mapping file. If the file looks like a NGS file (fasta, vcf, bam...) the mapping is extracted from a dictionary; Otherwise, it is interpreted as a mapping file ( See https://github.com/dpryan79/ChromosomeMappings )";
/** if file looks like a dictionary (fasta, vcf, dict...) use it , otherwise it's a mapping file */
public static ContigNameConverter fromPathOrOneDictionary(final Path file) {
	IOUtil.assertFileIsReadable(file);	
	final String filename=file.getFileName().toString();
	if(ReferenceSequenceFileFactory.FASTA_EXTENSIONS.stream().anyMatch(E->filename.endsWith(E)) || 
		StringUtils.endsWith(filename, 
				IOUtil.SAM_FILE_EXTENSION,
				BamFileIoUtils.BAM_FILE_EXTENSION,
				CramIO.CRAM_FILE_EXTENSION,
				IOUtil.DICT_FILE_EXTENSION,
				IOUtil.INTERVAL_LIST_FILE_EXTENSION,
				IOUtil.VCF_FILE_EXTENSION,
				IOUtil.COMPRESSED_VCF_FILE_EXTENSION,
				IOUtil.BCF_FILE_EXTENSION
				))
		{
		final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(file);
		return fromOneDictionary(dict);
		}
	else
		{
		return fromPath(file);
		}
	}


public static final String OPT_MAPPING_FILE_DESC="Chromosome mapping file. See https://github.com/dpryan79/ChromosomeMappings";

public static ContigNameConverter fromFile(final File mappingFile)
	{
	return fromPath(mappingFile==null?null:mappingFile.toPath());
	}


public static ContigNameConverter fromPath(final Path mappingFile)
	{
	IOUtil.assertFileIsReadable(mappingFile);
	final MapBasedContigNameConverter mapper=new MapBasedContigNameConverter();
	mapper.name=mappingFile.getFileName().toString();
	BufferedReader in=null;
	try {
		final CharSplitter tab= CharSplitter.TAB;
		in=IOUtils.openPathForBufferedReading(mappingFile);
		String line;
		while((line=in.readLine())!=null)
			{
			if(StringUtils.isBlank(line) || line.startsWith("#")) continue;
			final String tokens[]=tab.split(line);
			if(tokens.length!=2
					|| StringUtils.isBlank(tokens[0])
					|| StringUtils.isBlank(tokens[1])
					) {
				in.close();in=null;
				throw new IOException("Bad mapping line: \""+line+"\" in "+mappingFile);
				}
			tokens[0]=tokens[0].trim();
			tokens[1]=tokens[1].trim();
			if(mapper.map.containsKey(tokens[0]))
				{
				in.close();
				throw new IOException("Mapping defined twice for: \""+tokens[0]+"\"");
				}
			mapper.map.put(tokens[0], tokens[1]);
			}
		return mapper;
		}
	catch(final IOException err) {
		throw new RuntimeIOException(err);
		}
	finally {
		CloserUtil.close(in);
	}
}

public static ContigNameConverter fromMap(final Map<String,String> hash) {
	final MapBasedContigNameConverter mapper=new MapBasedContigNameConverter();
	mapper.name="fromMap";
	mapper.map.putAll(hash);
	return mapper;
	}


public static ContigNameConverter createConvertToUcsc()
	{
	final MapBasedContigNameConverter map = new MapBasedContigNameConverter();
	for(int i=1;i<=22;++i) map.put(String.valueOf(i),"chr"+i);
	map.put("X","chrX");map.put("Y","chrY");map.put("MT","chrM");
	return map;
	}


public static ContigNameConverter createConvertToEnsembl()
	{
	final MapBasedContigNameConverter map = new MapBasedContigNameConverter();
	for(int i=1;i<=22;++i) map.put("chr"+i,String.valueOf(i));
	map.put("chrX","X");map.put("chrY","Y");map.put("chrM","MT");
	return map;
	}

/*
private static class GRCh37Ucsc2Ensembl extends MapBasedContigNameConverter
	{
	GRCh37Ucsc2Ensembl() {
		for(int i=1;i<=22;++i) this.put("chr"+i, String.valueOf(i));
		this.put("chrX","X");this.put("chrY","Y");this.put("chrM","MT");
		}
	@Override
	public String getName() {return "GRCh37UcscToEnsembl";}
	}*/

/** creates a ContigNameConverter from two dictionaries, just using the common chromosome names */
public static ContigNameConverter fromDictionaries(
		final SAMSequenceDictionary dictIn,
		final SAMSequenceDictionary dictOut
		)
	{
	return new TwoDictionaries(dictIn,dictOut);
	}

/** creates a ContigNameConverter from one dictionary, try to guess best conversion*/
public static ContigNameConverter fromOneDictionary(
		final SAMSequenceDictionary dict
		)
	{
	return new OneDictionary(dict);
	}
/** creates a ContigNameConverter from a set of contig names, try to guess best conversion*/
public static ContigNameConverter fromContigSet(
		final Set<String> contigSet
		)
	{
	return new OneDictionary(contigSet);
	}

/** creates a ContigNameConverter from an IntervalTreeMap*/
public static ContigNameConverter fromIntervalTreeMap(
		final IntervalTreeMap<?> itm
		)
	{
	return fromContigSet(itm.keySet().stream().map(I->I.getContig()).collect(Collectors.toSet()));
	}


private static class TwoDictionaries extends MapBasedContigNameConverter
	{
	TwoDictionaries(final SAMSequenceDictionary dictIn,final SAMSequenceDictionary dictOut) {
		if(dictIn==null || dictIn.isEmpty())
			{
			throw new IllegalArgumentException("dictIn is null or empty");
			}
		if(dictOut==null || dictOut.isEmpty())
			{
			throw new IllegalArgumentException("dictOut is null or empty");
			}
		
		for(final SAMSequenceRecord ssr : dictIn.getSequences())
			{
			if(dictOut.getSequence(ssr.getSequenceName())!=null)
				{
				super.map.put(ssr.getSequenceName(), ssr.getSequenceName());
				continue;
				}
			
			if(ssr.getSequenceName().startsWith("chr") )
				{
				final String ctg = ssr.getSequenceName().substring(3);
				if(dictOut.getSequence(ctg)!=null)
					{
					super.map.put(ssr.getSequenceName(), ctg);
					continue;
					}
				if(ssr.getSequenceName().equals("chrMT") && dictOut.getSequence("M")!=null)
					{
					super.map.put("chrMT","M");
					continue;
					}
				if(ssr.getSequenceName().equals("chrM") && dictOut.getSequence("MT")!=null)
					{
					super.map.put("chrM","MT");
					continue;
					}
				continue;
				}
			if(!ssr.getSequenceName().startsWith("chr") )
				{
				final String ctg = "chr"+ssr.getSequenceName();
				if(dictOut.getSequence(ctg)!=null)
					{
					super.map.put(ssr.getSequenceName(), ctg);
					continue;
					}
				if(ssr.getSequenceName().equals("MT") && dictOut.getSequence("chrM")!=null)
					{
					super.map.put("MT","chrM");
					continue;
					}
				if(ssr.getSequenceName().equals("M") && dictOut.getSequence("chrMT")!=null)
					{
					super.map.put("M","chrMT");
					continue;
					}
				continue;
				}
			
			final String md5in = ssr.getMd5();
			if(!StringUtil.isBlank(md5in)) {
				final List<SAMSequenceRecord> ssrOUts = dictOut.getSequences().
						stream().
						filter(SSR->!StringUtil.isBlank(SSR.getMd5()) && md5in.equals(SSR.getMd5())).
						collect(Collectors.toList());
				if(ssrOUts.size()==1)
					{
					super.map.put(ssr.getSequenceName(), ssrOUts.get(0).getSequenceName());
					continue;
					}
				}
			}
		}
	@Override
	public String getName() {return "TwoDictionaries";}
	}


/** creates a ContigNameConverter from only one contig name. e.g: in LowResBam when we plot in one region. */
public static ContigNameConverter fromOneContig(final String contig)
	{
	if(StringUtil.isBlank(contig)) throw new IllegalArgumentException("empty contig");
	final Map<String,String> hash= new HashMap<>();
	hash.put(contig, contig);
	if(contig.startsWith("chr"))
		{
		hash.put(contig.substring(3), contig);
		}
	else if(!contig.startsWith("chr"))
		{
		hash.put("chr"+contig, contig);
		}
	if(contig.equals("MT")) hash.put("chrM", contig);
	if(contig.equals("chrM")) hash.put("MT", contig);
	return fromMap(hash);
	}


private static class OneDictionary extends ContigNameConverter
	{
	private final SAMSequenceDictionary dict;
	@SuppressWarnings("serial")
	private final Set<String> mitochrondrials = new HashSet<String>() {{{
		add("M");
		add("MT");
		add("chrM");
		add("chrMT");
	}}};
	OneDictionary( final SAMSequenceDictionary dict)
		{
		this.dict = dict;
		}
	
	OneDictionary( final Set<String> contigs)
		{
		this.dict = new SAMSequenceDictionary(
			contigs.stream().
				map(C->new SAMSequenceRecord(C, 9999)).
				collect(Collectors.toList())
			);
		}
	
	@Override
	protected String find(final String contig) {
		if(this.dict.getSequenceIndex(contig)!=-1) return contig;
		
		if(this.mitochrondrials.contains(contig))
			{
			final Optional<String> mitName = this.mitochrondrials.
				stream().
				filter(S->!S.equals(contig)).
				map(S->this.dict.getSequence(S)).
				filter(SSR->SSR!=null).
				map(SSR->SSR.getSequenceName()).
				findFirst();
			if(mitName.isPresent()) return mitName.get();
			}
		
		
		final String c2;
		if(contig.startsWith("chr"))
			{
			c2 = contig.substring(3);
			}
		else
			{
			c2 = "chr"+contig;
			}
		final SAMSequenceRecord ssr = this.dict.getSequence(c2);
		if(ssr!=null) return ssr.getSequenceName();
		
		// chrUn_gl000247 vs GL000247.1
		if(contig.startsWith("chrUn_gl"))
			{
			String c3=contig.substring(8);
			for(final SAMSequenceRecord ssr2 :this.dict.getSequences())
				{
				String c4 = ssr2.getSequenceName();
				if(!c4.startsWith("GL")) continue;
				c4=c4.substring(2);
				if(c4.equals(c3)) return ssr2.getSequenceName();
				if(!c4.startsWith(c3)) continue;
				c4=c4.substring(c3.length());
				if(c4.startsWith(".")) return ssr2.getSequenceName();
				}	
			}
		
		return null;
		}
	@Override
	public String getName() {
		if(SequenceDictionaryUtils.isGRCh37(this.dict)) return "GRCh37";
		if(SequenceDictionaryUtils.isGRCh38(this.dict)) return "GRCh38";
		return "OneDictionary";
		}
	}
/** add common known contig aliases to the dict.  eg chr1->1, 2 -> chr2
 * @return the dict
 */
public static SAMSequenceDictionary setDefaultAliases(final SAMSequenceDictionary dict ) {
	final Set<String> contigs = dict.getSequences().stream().map(C->C.getSequenceName()).collect(Collectors.toSet());
	final BiConsumer<String, String> alias = (S1,S2)->{
		if(contigs.contains(S2)) return;
		dict.addSequenceAlias(S1, S2);
		};
	final Predicate<String> isNumber = (S)->!S.isEmpty() && 
			S.chars().allMatch(ASCII->Character.isDigit(ASCII));
		
	for(final String C: contigs)
		{
		if(C.startsWith("chr"))
			{
			if(C.equals("chrX")) {
				alias.accept(C, "X");
				}
			else if(C.equals("chrY")) {
				alias.accept(C, "Y");
				}
			else if(C.equals("chrM") || C.equals("chrMT")) {
				alias.accept(C, "MT");
				alias.accept(C, "M");
				}
			else if(isNumber.test(C.substring(3))) {
				alias.accept(C,C.substring(3));
				}
			}
		else
			{
			if(C.equals("X")) {
				alias.accept(C,"chrX");
				}
			else if(C.equals("Y")) {
				alias.accept(C,"chrY");
				}
			else if(C.equals("M") || C.equals("MT")) {
				alias.accept(C,"chrM");
				alias.accept(C,"chrMT");
				}
			else if(isNumber.test(C)) {
				alias.accept(C,"chr"+C);
				}
			}
		}
	return dict;
	}
}
