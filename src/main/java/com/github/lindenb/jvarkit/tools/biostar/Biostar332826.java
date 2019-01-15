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
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.vcf.VCFConstants;
/**
 BEGIN_DOC
 
 ## Example
 
 ```
 $ wget -O - -q "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" |\
 	gunzip -c |\
 	java -jar dist/biostar332826 --ids ids.txt > out.vcf 
 ```
 
 END_DOC
 */
@Program(name="biostar332826",
description="Fast Extraction of Variants from a list of IDs",
keywords= {"vcf","rs"},
biostars=332826
)
public class Biostar332826 extends Launcher {
private static final Logger LOG = Logger.build(Biostar332826.class).make();

@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
private File outputFile = null;
@Parameter(names={"-r","-i","--ids"},description="A list of identifiers, one per line")
private File rsFile = null;
@Parameter(names={"-R","-I"},description="A semicolon/comma/space separated list of identifiers")
private String rsStr = "";
@Parameter(names={"-d","--delete"},description="When found , remove the ID from the list of identifiers. Should be faster but don't use it if two variants have the same ID.")
private boolean removeIfFound=false;
@Parameter(names={"-v","--inverse"},description="Inverse: don't print the variants containing the IDS.")
private boolean inverseSelection =false;


@Override
public int doWork(final List<String> args) {
	final Set<String> rsSet = new HashSet<>();
	BufferedReader br = null;
	PrintWriter pw = null;
	try {
		if(this.rsFile!=null)
			{
			rsSet.addAll(IOUtil.slurpLines(this.rsFile));
			}	
		for(final String str:this.rsStr.split("[ ;,]"))
			{
			if(StringUtil.isBlank(str)) continue;
			rsSet.add(str);
			}
		LOG.info("rs list size: "+rsSet.size());
		rsSet.remove("");
		if(rsSet.isEmpty()) LOG.warn("NO IDENTIFIER WAS SPECIFIED");
		br = super.openBufferedReader(oneFileOrNull(args));
		pw = super.openFileOrStdoutAsPrintWriter(this.outputFile);
		String line;
		while((line=br.readLine())!=null) {
			if(!line.startsWith("#")) {
				LOG.error("VCF header line doesn't start with '#' "+line);
				return -1;
				}
			pw.println(line);
			if(line.startsWith("#CHROM")) break;
			}
		
		final CharSplitter tab = CharSplitter.TAB;
		
		if(this.inverseSelection)
			{
			while((line=br.readLine())!=null) {
				final List<CharSequence> tokens = tab.splitAsCharSequenceList(line, 4);
				if(tokens.size() != 4) {
					LOG.error("expected at least four tokens in "+line);
					return -1;
					}
				final String id = tokens.get(2).toString();
				
				if(!(this.removeIfFound && !id.equals(VCFConstants.EMPTY_ID_FIELD) ?
						rsSet.remove(id):rsSet.contains(id)))
					{
					pw.println(line);
					}
				}
			IOUtils.copyTo(br, pw);
			}
		else
			{
			while((line=br.readLine())!=null && !rsSet.isEmpty()) {
				final List<CharSequence> tokens = tab.splitAsCharSequenceList(line, 4);
				if(tokens.size() !=4) {
					LOG.error("expected at least four tokens in "+line);
					return -1;
					}
				final String id = tokens.get(2).toString();
				if(id.isEmpty()) continue;
				if(this.removeIfFound && !id.equals(VCFConstants.EMPTY_ID_FIELD) ? rsSet.remove(id):rsSet.contains(id))
					{
					pw.println(line);
					}	
				}
			}
		pw.flush();
		pw.close();pw = null;
		br.close();br=null;
		return 0;
		}
	catch(final Exception err) {
		LOG.error(err);
		return -1;
		}
	finally
		{
		CloserUtil.close(br);
		CloserUtil.close(pw);
		}
	}	
public static void main(final String[] args) throws IOException
	{
	new Biostar332826().instanceMainWithExit(args);
	}
}
