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
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

/** fasta reference file profide. */
public abstract class ReferenceFileSupplier implements Supplier<File> {
	private static final Logger LOG = Logger.build(ReferenceFileSupplier.class).make();
	/** java propery */
	private static final String JAVA_PROP = "jvarkit.fasta.reference";
	/** linux env variabl */
	private static final String ENV_KEY = "FASTA_REFERENCE";
	/** catalog file */
	private static final String CATALOG_FILE ="fasta-ref.properties";
	
	
	private static final String CatalogKeyPatStr = "[A-Za-z][A-Za-z0-9_\\\\-]*";
	private static final Pattern CatalogKeyPattern = Pattern.compile(CatalogKeyPatStr);
	
	/** description for jcommander args */
	public static final String OPT_DESCRIPTION =
			"The parameter is the path to an Indexed fasta Reference file. "+
			"This fasta file must be indexed with samtools faidx and with picard CreateSequenceDictionary. " +
			"The parameter can also be a 'key' (matching the regular expression `" + CatalogKeyPatStr +"`) in a catalog file. " +
			"A 'catalog' file is a java property file ( https://docs.oracle.com/javase/tutorial/essential/environment/properties.html ) where the values are the path to the fasta file. " +
			" Catalogs are searched in that order : `${PWD}/"+CATALOG_FILE +"`, `${HOME}/."+CATALOG_FILE +"`, `/etc/jvarkit/" + CATALOG_FILE +"`. " +
			" If the key or the path are not defined by the user, they will be searched in that order 1) the java property -D"+JAVA_PROP+"=pathTofastaOrCatalogKey ." +
			" 2) the linux environement variable $"+ENV_KEY+"=pathTofastaOrCatalogKey 3) The catalogs."
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
						"Cannot find key \""+path+"\" in the Fasta Reference Catalogs / property files ( "+ 
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
		public synchronized File get() {
			if(this.searched) {
				return this.fasta;
				}
			this.searched = true;
			this.fasta = parseVal(System.getProperty(JAVA_PROP,null));
			if(this.fasta!=null) return fasta;
			this.fasta = parseVal(System.getenv(ENV_KEY));
			if(this.fasta!=null) return fasta;
			LOG.warn("cannot find reference fasta file using key/path = \""+this.fasta+"\".\n" +
					"\t${"+ENV_KEY+"}="+System.getenv(ENV_KEY)+"\n" +
					"\t-D"+JAVA_PROP+"="+System.getProperty(JAVA_PROP,null)+"\n"+
					getCatalogFiles().stream().map(F->"\tCatalog file: "+ F.getPath()+" exists:"+F.exists()).collect(Collectors.joining("\n")) + "\n"+
					"Check the application parameters please.");
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
	
	/** list all potential catalog files, even if they don't exist */
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
	
	/** search catalogs using a samsequencedictionary. Return NULL if not found or the dict is null*/
	@SuppressWarnings("deprecation")
	public static File searchCatalogs(final SAMSequenceDictionary dict) {
		if(dict==null) return null;
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
					final SAMSequenceDictionary d2 = SAMSequenceDictionaryExtractor.extractDictionary(F);
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
		return getCatalogFiles().
				stream().
				map(C->searchCatalog(C, key)).
				filter(F->F!=null).
				findFirst().
				orElse(null);
		}
	
	/** return the path to the fasta file. May be null */
	@Override
	public abstract File get();
	
	@Override
	public String toString() {
		return "Fasta Reference File: "+get();
		}
	}
