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
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.table.Table;
import com.github.lindenb.jvarkit.table.TableFactory;
import com.github.lindenb.jvarkit.table.TextExporter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;

@Program(name="prettytable",
description="Pretty tabular data",
keywords={"util","table","tsv","tabular","data"},
generate_doc=false
)
public class PrettyTable extends Launcher{
	private static final Logger LOG = Logger.build(PrettyTable.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-d","--delim"},description="char delimiter.")
	private String delimiter="\\t";
	@Parameter(names={"-t","--transpose"},description="transpose matrix")
	private boolean transpose = false;
	@Parameter(names={"-n","--no-header"},description="The is no header")
	private boolean noHeader = false;
	@Parameter(names={"-u","--unicode"},description="Allow unicode characters")
	private boolean withUnicode = false;

	
@Override
public int doWork(final List<String> args) {
	final List<List<String>> rows=new ArrayList<>();
	PrintWriter out=null;
	BufferedReader br=null;
	final CharSplitter delim = CharSplitter.of(this.delimiter);
	try {
		int n_columns=-1;
		br= super.openBufferedReader(oneFileOrNull(args));
		String line;
		while((line=br.readLine())!=null)
			{
			if(StringUtils.isBlank(line)) continue;
			final List<String> row = delim.splitAsStringList(line);
			if(rows.isEmpty())
				{
				n_columns = row.size();
				if(this.noHeader)
					{
					final List<String> header = new ArrayList<>(n_columns);
					while(header.size()< n_columns)
						{
						header.add("$"+(1+header.size()));
						}
					rows.add(header);
					}
				}
			else if(row.size()<n_columns)
				{
				while(row.size()<n_columns) row.add("");
				}
			else if(row.size()>n_columns)
				{
				n_columns = row.size();
				for(int y=0;y< rows.size();++y)
					{
					while(rows.get(y).size()<n_columns) rows.get(y).add("");
					}
				}
			else
				{
				//ok
				}
			rows.add(row);
			}
		br.close();
		br=null;
		
		if(rows.isEmpty()) {
			LOG.error("empty input");
			return -1;
			}
		
		if(this.transpose)
			{
			final List<List<String>> rows2 = new ArrayList<>(n_columns);
			for(int x=0;x< n_columns;++x)
				{
				final List<String> row = new ArrayList<>(rows.size());
				for(int y=0;y< rows.size();++y)
					{
					row.add(rows.get(y).get(x));
					}
				rows2.add(row);
				}
			rows.clear();
			rows.addAll(rows2);
			}
		
		final TableFactory tableFactory = new TableFactory();
		final Table table = tableFactory.createTable(rows.get(0));
		for(int y=1;y< rows.size();++y)
			{	
			table.addRow(rows.get(y));
			}
		
		out= super.openFileOrStdoutAsPrintWriter(this.outputFile);
		TextExporter exporter=new TextExporter();
		exporter.setAllowUnicode(this.withUnicode);
		
		exporter.print(table, out);
		out.flush();
		out.close();
		out = null;
		return 0;
		}
	catch(final Exception err) {
		LOG.error(err);
		return -1;
		}
	finally
		{
		CloserUtil.close(out);
		}
	}
	
public static void main(final String[] args) {
	new PrettyTable().instanceMainWithExit(args);
	}
}
