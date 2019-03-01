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
package com.github.lindenb.jvarkit.dict;

import java.io.IOException;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.Properties;
import java.util.function.Predicate;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

/** 
 * Registry used to find REFerence sequences
 *
 */
public interface ReferenceRegistry {
	static final Logger LOG = Logger.build(ReferenceRegistry.class).make();
	/** catalog file */
	static final String CATALOG_FILE ="fasta-ref.properties";

	/** description for jcommander args */
	public static final String REGISTRY_DESCRIPTION =
			"A 'catalog' file is a java property file ( https://docs.oracle.com/javase/tutorial/essential/environment/properties.html ) where the values are the path to the fasta file. " +
			"Catalogs are searched in that order : `${PWD}/"+CATALOG_FILE +"`, `${HOME}/."+CATALOG_FILE +"`, `/etc/jvarkit/" + CATALOG_FILE +"`. " +
			"A fasta file must be indexed with samtools faidx and with picard CreateSequenceDictionary."
			;
/** find Path to reference sequence using a dictionary */
public Optional<Path> getReferenceByDictionary(final SAMSequenceDictionary dict);
/** find Path to reference sequence using a name. e.g "hg19" */
public Optional<Path> getReferenceByName(final String name);

/** find Dict of reference sequence using a name. e.g "hg19" */
public default Optional<SAMSequenceDictionary> getDictionaryByName(final String name) {
	Optional<Path> path = this.getReferenceByName(name);
	if(!path.isPresent()) return Optional.empty();
	return Optional.of(SequenceDictionaryUtils.extractRequired(path.get()));
}


/** find Path to reference sequence using the dictionary of a dict, vcf or a bam file*/
public default Optional<Path> getReferenceByPath(final Path path) {
	try {
		final SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(path);
		if(dict==null) return Optional.empty();
		return getReferenceByDictionary(dict);
	} catch(final Throwable err) {
		return Optional.empty();
		}
	}


public static	ReferenceRegistry  getDefault() {
	return new RegistryImpl();
	}

static class RegistryImpl implements ReferenceRegistry
	{
	private static class Entry
		{
		final String name;
		final  Path ref;
		Entry(String name,Path ref) {
			this.name = name;
			this.ref = ref;
			}
		}
	
	private boolean isFile(final Path path) {
		return path!=null && Files.exists(path) && Files.isRegularFile(path)&& Files.isReadable(path);
		}
	
	private Path asValidRef(final Path fasta) {
		if(!isFile(fasta)) {
			LOG.warn("No file \""+fasta+"\"");
			return null;
			}
		
		final Path fai = ReferenceSequenceFileFactory.getFastaIndexFileName(fasta);
		if(!isFile(fai)) {
			LOG.warn("No .fai file for "+fasta+". Use `samtools faidx` to create this file.");
			return null;
		}

		final Path dict=ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(fasta);
		if(!isFile(dict)) {
			LOG.warn("No .dict file for "+fasta+". Use `picard CreateSequenceDictionary` to create this file.");
			return null;
		}
		
		return fasta;
		}
	/** search catalog file , return all entry with a valid fasta file */
	private List<Entry> listEntries(final Path propFile,final Predicate<Entry> predicate) {
		if(!isFile(propFile)) return Collections.emptyList();

		
		final Properties props = new Properties();
		try (final Reader r = Files.newBufferedReader(propFile)) {
			props.load(r);
			}
		catch(final IOException err)
			{
			return Collections.emptyList();
			}
		final List<Entry> L = new ArrayList<>(props.size());
		for(final Object k: props.keySet()) {
			final String key = k.toString();
			if(StringUtils.isBlank(key)) return null;
			final String value = props.getProperty(key);
			if(StringUtils.isBlank(value)) return null;
			final Entry entry = new Entry(key,asValidRef(Paths.get(value)));
			if(entry.ref==null) continue;
			if(predicate!=null && !predicate.test(entry)) continue;
			L.add(entry);
			}
		return L;
		}
	
	private List<Path> getRegistries() {
		final List<Path> list = new ArrayList<>();
		/* user cwd */
		final Path catFile1 = Paths.get(System.getProperty("user.dir",null),CATALOG_FILE);
		if(isFile(catFile1)) list.add(catFile1);
		
		/* user home */
		final  Path catFile2 = Paths.get(System.getProperty("user.home",null),"." + CATALOG_FILE);
		if(isFile(catFile2)) list.add(catFile2);
		
		/* etc */
		final  Path catFile3 = Paths.get("/etc/jvarkit" , CATALOG_FILE);
		if(isFile(catFile3)) list.add(catFile3);
		
		if(list.isEmpty()) {
			LOG.warn("none of the following reference registry was found: "+
					catFile1+" "+catFile2+" "+catFile3+"\n"+
					REGISTRY_DESCRIPTION
					);
			}
		
		return list;
		}
	
	@Override
	public Optional<Path> getReferenceByDictionary(final SAMSequenceDictionary dict) {
		final Predicate<Entry> pred= P->SequenceUtil.areSequenceDictionariesEqual(dict, SequenceDictionaryUtils.extractRequired(P.ref));
		return getRegistries().
			stream().
			flatMap(F->listEntries(F,pred).stream()).
			map(E->E.ref).
			findFirst();
		}

	@Override
	public Optional<Path> getReferenceByName(final String name) {
		final Predicate<Entry> pred= P->P.name.equals(name);
		return getRegistries().
			stream().
			flatMap(F->listEntries(F,pred).stream()).
			map(E->E.ref).
			findFirst();
		}
	
	
	}



}
