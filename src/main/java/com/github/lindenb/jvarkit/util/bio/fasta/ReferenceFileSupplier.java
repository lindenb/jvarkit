/*
The MIT License (MIT)

Copyright (c) 2018 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.util.bio.fasta;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;
import java.util.function.Supplier;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.beust.jcommander.IStringConverter;
import com.beust.jcommander.ParameterException;
import com.github.lindenb.jvarkit.lang.JvarkitException;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

/** fasta reference file profide. */
public abstract class ReferenceFileSupplier implements Supplier<File> {
	/** java propery */
	private static final String JAVA_PROP = "jvarkit.fasta.reference";
	/** linux env variabl */
	private static final String ENV_KEY = "FASTA_REFERENCE";
	/** catalog file */
	private static final String CATALOG_FILE ="fasta-ref.properties";
	
	private static final Pattern CatalogKeyPattern = Pattern.compile("[A-Za-z][A-Za-z0-9_\\\\-]*");
	public static final String OPT_DESCRIPTION =
			"Path to an Indexed fasta Reference file. "+
			"This fasta file must be indexed with samtools faidx and with picard CreateSequenceDictionary. " +
			"The parameter can also be a catalog key matching " + CatalogKeyPattern.pattern() +" . " +
			"A catalog is a java property file ( https://docs.oracle.com/javase/tutorial/essential/environment/properties.html ) where the values are the path to the fasta file. " +
			" Catalog are searched in that order : $‘PWD}/"+CATALOG_FILE +", ${HOME}/."+CATALOG_FILE +", /etc/jvarkit/" + CATALOG_FILE +
			"If empty the key or the path will be searched in that order 1) the java property -D"+JAVA_PROP+"=pathTofastaOrCatalogKey ." +
			" 2) the linux environement variable $"+ENV_KEY+"=pathTofastaOrCatalogKey 3) the catalogs."
			;

	
public static class StringConverter
	implements IStringConverter<ReferenceFileSupplier>
	{
	@Override
	public ReferenceFileSupplier convert(final String path) {
		if(!StringUtil.isBlank(path))
			{
			final File f;
			if(CatalogKeyPattern.matcher(path).matches()) {
				f = ReferenceFileSupplier.searchCatalogs(path);
				if(f==null) throw new ParameterException(
						"Cannot find key "+path+" in the Fasta Reference Catalogs ( "+ 
						getCatalogFiles().stream().map(F->F.getPath()).collect(Collectors.joining(" ")) +
						"). You might provide the path to a fasta reference instead of a catalog key"
						);
				}
			else
				{
				f = new File(path);
				}
			
			IOUtil.assertFileIsReadable(f);
			return new FileSupplier(f);
			}
		else
			{
			return getDefaultReferenceFileSupplier();
			}
		}
	}
	
	private static class FileSupplier extends ReferenceFileSupplier
		{
		final File fastaFile;
		
		private FileSupplier(final File fastaFile) {
			this.fastaFile = fastaFile;
			}
		
		@Override
		public File get() {
			return this.fastaFile;
			}
		
		}
	
	public static ReferenceFileSupplier getDefaultReferenceFileSupplier() {
		return new DefaultSupplier();
		}
	
	private static class DefaultSupplier extends ReferenceFileSupplier
		{
		private boolean searched = false;
		private File fasta=null;
		
		private File parseVal(final String s) {
			if(StringUtil.isBlank(s)) return null;
			if(CatalogKeyPattern.matcher(s).matches())
				{
				return searchCatalogs(s);
				}
			final File fasta = new File(s);
			if(fasta.exists() && fasta.isFile()) return fasta;
			return null;
			}
		
		@Override
		public File get() {
			if(this.searched) {
				return this.fasta;
				}
			this.searched = true;
			this.fasta = parseVal(System.getProperty(JAVA_PROP,null));
			if(this.fasta!=null) return fasta;
			this.fasta = parseVal(System.getenv(ENV_KEY));
			if(this.fasta!=null) return fasta;
			return null;
			}
		}
	
	/** get required fasta file or throws an exception */
	public File getRequired() {
		final File f = get();
		if(f!=null) return f;
		throw new JvarkitException.ReferenceMissing("Reference file is missing. ");
		}

	/** search catalog file for the given key */
	private static File searchCatalog(final File propFile,final String key) {
		if(StringUtil.isBlank(key) || propFile==null || !propFile.exists() || !propFile.isFile() || !propFile.canRead()) return null;
		final Properties props = new Properties();
		try (final FileReader r = new FileReader(propFile)) {
			props.load(r);
			}
		catch(final IOException err)
			{
			return null;
			}
		if(!props.containsKey(key)) return null;
		final String v = props.getProperty(key);
		if(StringUtil.isBlank(v)) return null;
		final File fasta = new File(v);
		if(fasta.exists() && fasta.canRead() && fasta.isFile()) return fasta;
		return null;
		}
	
	private static List<File> getCatalogFiles() {
		final List<File> list = new ArrayList<>();
		/* user cwd */
		File catFile = new File(System.getProperty("user.dir",null),CATALOG_FILE);
		list.add(catFile);
		
		/* user home */
		catFile = new File(System.getProperty("user.home",null),"." + CATALOG_FILE);
		list.add(catFile);
		
		/* etc */
		catFile = new File("/etc/jvarkit/" + CATALOG_FILE);
		list.add(catFile);
		
		return list;
		}
	
	/** search catalogs using a samsequencedictionary. Return NULL if not found */
	public static File searchCatalogs(final SAMSequenceDictionary dict) {
		for(final File catFile: getCatalogFiles()) {
			if(!catFile.exists() || !catFile.isFile() || !catFile.canRead()) continue;
			
			final Properties props = new Properties();
			try (final FileReader r = new FileReader(catFile)) {
				props.load(r);
				}
			catch(final IOException err)
				{
				continue;
				}
			final File fasta = props.values().stream().
				filter(O->O!=null).
				map(O->O.toString()).
				map(S->new File(S)).
				filter(F->F.exists() && F.isFile() && F.canRead()).
				filter(F->{
					final SAMSequenceDictionary d2=SAMSequenceDictionaryExtractor.extractDictionary(F);
					if(d2==null || !SequenceUtil.areSequenceDictionariesEqual(d2, dict)) return false;
					return true;
					}).
				filter(F->F!=null).
				findFirst().orElse(null);
			if(fasta!=null) return fasta;
			}
		return null;
		}
	
	public static File searchCatalogs(final String key) {
		for(final File f: getCatalogFiles()) {
			File fasta =  searchCatalog(f,key);
			if(fasta!=null) return fasta;
			}
		return null;
		}
	
	@Override
	public String toString() {
		return "Fasta Reference File: "+get();
		}
	}
