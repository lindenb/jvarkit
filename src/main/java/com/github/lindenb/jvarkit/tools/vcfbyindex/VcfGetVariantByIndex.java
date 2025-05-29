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
package com.github.lindenb.jvarkit.tools.vcfbyindex;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.OptionalLong;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;
import com.github.lindenb.jvarkit.variant.vcf.VCFByIndex;

import htsjdk.tribble.util.ParsingUtils;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
/**
BEGIN_DOC

## Example

```bash
# get random indexes
$ gunzip -c input.vcf.gz | grep -v "#" | awk '{print NR;}' | shuf | head -n 15 > index.list
# get those 15 variants
java -jar dist/jvarkit.jar vcfgetvariantbyindex --force --index-file index.list input.vcf.gz > output.vcf

```

END_DOC
*/
@Program(name="vcfbyindex",
	description="Access a Plain or BGZF-compressed VCF file by index",
	keywords={"vcf"},
	modificationDate = "20250529",
	jvarkit_amalgamion = true,
	menu="VCF Manipulation"
	)
public class VcfGetVariantByIndex extends Launcher
	{
	private static Logger LOG=Logger.of(VcfGetVariantByIndex.class);
	
	@Parameter(names="-o",description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;
	@Parameter(names={"-I","--index-file"},description=" (file) list of 1-based indexes for query")
	private Path fileListOfIndexes=null;
	@Parameter(names={"-i","--index"},description="Comma separated of 1-based index for query")
	private String indexStr="";
	@Parameter(names={"-f","--force"},description="Force the creation of the index")
	private boolean force_index_creation=false;
	@ParametersDelegate
	private WritingVariantsDelegate writingVariantsDelegate= new WritingVariantsDelegate();

	

	private static OptionalLong parseIndex(final VCFByIndex index, final String line) {
		long ith;
		try {
			ith=Long.parseLong(line);
			} 
		catch (NumberFormatException e) {
			LOG.error("Bad index in "+line+" ignoring");
			return OptionalLong.empty();
			}
		ith--;//0-based index
		if(ith<0 || ith>=index.size())
			{
			LOG.warn("Index out of bound in "+line+" ignoring");
			return OptionalLong.empty();
			}
		return OptionalLong.of(ith);
		}
	
	public int doWork(final List<String> args) {
		try {
			if(args.size()!=1 && args.size()!=2) {
				LOG.error("illegal number of arguments. expect : VCF or VCF IDX");
				return -1;
				}
			final String vcfSource = args.get(0);
			final String indexSource = args.size()==1?
						ParsingUtils.appendToPath(vcfSource, VCFByIndex.INDEX_SUFFIX):
						args.get(1)
						;
			if(!VCFByIndex.indexExists(indexSource)) {
				if( this.force_index_creation) {
					if(IOUtils.isRemoteURI(indexSource)) {
						LOG.error("cannot create index "+indexSource+" for remote URL");
						return -1;
						}
					LOG.info("writing index "+indexSource);
					VCFByIndex.buildIndex(vcfSource, Paths.get(indexSource));
					}
				else
					{
					LOG.info("index not found "+indexSource+". use --force to create one.");
					return -1;
					}
				}
			
			try(VCFByIndex indexFile = new VCFByIndex(vcfSource, indexSource)) {
				try(VariantContextWriter w = this.writingVariantsDelegate.dictionary(indexFile.getHeader()).open(outputFile)) {
					final VCFHeader h2=new VCFHeader(indexFile.getHeader());
					JVarkitVersion.getInstance().addMetaData(this, h2);
					w.writeHeader(h2);
					if(fileListOfIndexes!=null) {
						try(BufferedReader r=IOUtils.openPathForBufferedReading(fileListOfIndexes)) {
							String line;
							while((line=r.readLine())!=null)
								{
								if(StringUtils.isBlank(line) || line.startsWith("#")) continue;
								final OptionalLong ith = parseIndex(indexFile,line);
								if(ith.isEmpty()) continue;
								w.add(indexFile.get(ith.getAsLong()));
								}
							}
						}
					if(!StringUtils.isBlank(this.indexStr)) {
						for(String stri : CharSplitter.COMMA.split(this.indexStr)) {
							if(StringUtils.isBlank(stri)) continue;
							final OptionalLong ith = parseIndex(indexFile,stri);
							if(ith.isEmpty()) continue;
							w.add(indexFile.get(ith.getAsLong()));
							}
						}
					}
				}
			return 0;
			} 
		catch(final Throwable e)
			{
			LOG.error(e);
			return -1;
			}
		}
	

	
	public static void main(final String[] args) throws IOException
		{
		new VcfGetVariantByIndex().instanceMainWithExit(args);
		}
	}
