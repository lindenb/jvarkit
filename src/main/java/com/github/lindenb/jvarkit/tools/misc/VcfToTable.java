/*
The MIT License (MIT)

Copyright (c) 2017 Pierre Lindenbaum

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

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

/**

 */
@Program(name="vcf2table",
		description="convert a vcf to a table, to ease display in the terminal",
		keywords={"vcf","table","visualization"})
public class VcfToTable extends Launcher {
	private static final Logger LOG = Logger.build(VcfToTable.class).make();

	private static class Column
		{
		final String label;
		int maxLength=0;
		Column(final String label) {
			this.label=label;
			this.maxLength=label.length();
			}
		}

	
	private static class Table
		{
		private String caption="";
		private final List<Column> columns;
		private final List<List<Object>> rows = new ArrayList<>();
		public Table(final String...header) {
			this(Arrays.asList(header));
			}
		public Table(final List<String> labels)
			{
			this.columns=labels.stream().
					map(L->new Column(L)).
					collect(Collectors.toList());
			}
		public Table setCaption(final String t) {
			this.caption=t;
			return this;
		}
		public void addRow(final Object...items)
			{	
			addList(Arrays.asList(items));
			}
		public void addList(final List<Object> row)
			{
			rows.add(row);
			for(int i=0;i< row.size()&& i< this.columns.size() ;++i)
				{
				final Object o = row.get(i);
				if(o==null) continue;
				final String str=o.toString();
				final int len = str.length();
				this.columns.get(i).maxLength = Math.max(this.columns.get(i).maxLength,len);
				}
			}
		private int tableWidth() {
			return this.columns.stream().mapToInt(C->C.maxLength +2 ).sum() 
					+ 1
					+ this.columns.size();
			}
	
		
		public void print(final PrintStream out) {
			if(this.rows.isEmpty()) return;
			final int tw = this.tableWidth();
			StringBuilder hr= new StringBuilder();
			
			out.println(this.caption);
			hr.append('+');
			for(int i=0;i< this.columns.size();i++)
				{
				hr.append("-");
				for(int j=0;j< this.columns.get(i).maxLength;j++)hr.append('-');
				hr.append("-+");
				}
			out.println(hr.toString());
			
			out.print('|');
			for(int i=0;i< this.columns.size();i++)
				{
				out.print(" ");
				out.print(this.columns.get(i).label);
				for(int j=this.columns.get(i).label.length();j< this.columns.get(i).maxLength;j++) out.print(' ');
				out.print(" |");
				}
			out.println();
			out.println(hr.toString());
			
			for(int y=0;y< this.rows.size();++y) {
				final List<Object> row= this.rows.get(y);
				out.print("|");
				for(int i=0;i< this.columns.size();i++)
					{
					String str= row.get(i)==null?"":row.get(i).toString();
					out.print(" ");
					out.print(str);
					for(int j=str.length();j< this.columns.get(i).maxLength;j++) out.print(' ');
					out.print(" |");
					}
				out.println();
				
			}
			out.println(hr.toString());
					
			}
		}
	
	public static class TerminalViewer
		implements VariantContextWriter
		{
		@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
		private File outputFile = null;
		@Parameter(names={"-H"},description="Print Header")
		private boolean printHeader=false;
		private int countVariants=0;
		
		private PrintStream out= System.out;
		private VCFHeader header=null;
		TerminalViewer() {
			}
		@Override
		public void writeHeader(VCFHeader header) {
			this.header = header;
			if(outputFile!=null) {
				try {
					this.out = new PrintStream(IOUtils.openFileForWriting(this.outputFile));
				} catch (final IOException e) {
					throw new RuntimeIOException(e);
				}
				
				}
			if(printHeader)
				{
				Table t=new Table("ID","Description").setCaption("FILTERS");
				header.getFilterLines().forEach(
					L->t.addRow(L.getID(),L.getDescription())
					);
				t.print(out);
				out.println();
				
				}
			}
		@Override
		public void add(final VariantContext vc) {
			if(out==null) return;
			++countVariants;
			out.println(">>"+vc.getContig()+"/"+vc.getStart()+"/"+vc.getReference().getDisplayString()+" "+countVariants);
			
			Table t=new Table("Key","Value").setCaption("Variant");
			t.addRow("chrom",vc.getContig());
			t.addRow("start",vc.getStart());
			t.addRow("end",vc.getEnd());
			t.addRow("ID",vc.hasID()?vc.getID():".");
			t.addRow("REF",vc.getReference().getDisplayString());
			
			
			t.print(out);
			
			 t=new Table("Idx","REF","Allele").setCaption("Alleles");
			 for(final Allele a: vc.getAlleles())
			 	{
				t.addRow(vc.getAlleleIndex(a),a.isReference()?"*":"",a.getDisplayString());
			 	}
			t.print(out);

			 t=new Table("Filter").setCaption("FILTERS");
			 for(final String f:vc.getFilters())
			 	{
				t.addRow(f);
			 	}
			t.print(out);

			
			out.println("<<"+countVariants);
			out.println();
			}
		@Override
		public boolean checkError() {
			if(out==null) return true;
			return out.checkError();
			}
		@Override
		public void close() {
			if(out==null) return;
			out.flush();
			out.close();
			out=null;
			}
		}
	
	
	@ParametersDelegate
	private TerminalViewer  viewer = new TerminalViewer();

	@Override
	public int doWork(List<String> args) {
		VcfIterator in = null;
		
		
		try {
			in = super.openVcfIterator(oneFileOrNull(args));
			VCFUtils.copyHeaderAndVariantsTo(in, viewer);
			viewer.out.close();viewer.out=null;
			in.close();in=null;
			return 0;
		} catch (Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			CloserUtil.close(viewer);
			}
		}
	public static void main(String[] args) {
		new VcfToTable().instanceMainWithExit(args);
	}
}
