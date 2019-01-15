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
* 2014 creation

*/
package com.github.lindenb.jvarkit.util.vcf.swing;

import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.regex.Pattern;

import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;

import com.github.lindenb.jvarkit.tools.vcfviewgui.VcfViewGui;

public class InfoTreeModel extends DefaultTreeModel
	{
	private static Log LOG=Log.getInstance(VcfViewGui.class);
	private Pattern pipeVep=Pattern.compile("[\\|]");
	private Pattern pipeSnpEff=Pattern.compile("[\\|\\(\\)]");

	private static final long serialVersionUID = 1L;
	public InfoTreeModel()
		{
		super(new DefaultMutableTreeNode("INFO",true));
		}
	
	private DefaultMutableTreeNode getTreeNodeRoot()
		{
		return DefaultMutableTreeNode.class.cast(getRoot());
		}	
	
	public void setContext(VariantContext ctx,VCFHeader header)
		{
		getTreeNodeRoot().removeAllChildren();
		if(ctx!=null)
			{
			for(String key:ctx.getAttributes().keySet())
				{
				Object v=ctx.getAttribute(key);
				Object o[];
				if(v==null)
					{
					o=new Object[]{null};
					}
				else if(v instanceof java.util.Collection)
					{
					o=((java.util.Collection<?>)v).toArray();
					}
				else if(v.getClass().isArray())
					{
					o=(Object[])v;
					}
				else
					{
					o=new Object[]{v};
					}
				if(o.length==0)
					{
					
					}
				else if(o.length==1)
					{
					if(key.equals("CSQ") || key.equals("EFF"))
						{
						List<String> columns= key.equals("CSQ")?getCSQCols(header):getEFFCols(header);
						DefaultMutableTreeNode n1=new DefaultMutableTreeNode(
								"<html><b>"+key+"</html>"
								,true);
						getTreeNodeRoot().add(n1);
						String tokens[]=(key.equals("CSQ")?pipeVep:pipeSnpEff).split(String.valueOf(o));
						for(int i=0;i< tokens.length && i< columns.size();++i)
							{
							DefaultMutableTreeNode n2=new DefaultMutableTreeNode(
									"<html><b>"+columns.get(i)+"</b>:"+tokens[i]+"</html>"
									,false);
							n1.add(n2);
							}
						}
					else
						{
						DefaultMutableTreeNode n=new DefaultMutableTreeNode(
								"<html><b>"+key+"</b>:"+o[0]+"</html>"
								,false);
						getTreeNodeRoot().add(n);
						}
					}
				else 
					{
					if(key.equals("CSQ") || key.equals("EFF"))
						{
					
						List<String> columns= key.equals("CSQ")?getCSQCols(header):getEFFCols(header);
						DefaultMutableTreeNode n1=new DefaultMutableTreeNode(
								"<html><b>"+key+"</html>"
								,true);
						getTreeNodeRoot().add(n1);
						int index=0;
						for(Object v2:o)
							{
							DefaultMutableTreeNode n2=new DefaultMutableTreeNode(
									String.valueOf(++index),
									true
									);
							n1.add(n2);
							String tokens[]= (key.equals("CSQ")?pipeVep:pipeSnpEff).split(String.valueOf(v2));
							for(int i=0;i< tokens.length && i< columns.size();++i)
								{
								DefaultMutableTreeNode n3=new DefaultMutableTreeNode(
										"<html><b>"+columns.get(i)+"</b>:"+tokens[i]+"</html>"
										,false);
								n2.add(n3);
								}
							}
						}
					else
						{
						DefaultMutableTreeNode n1=new DefaultMutableTreeNode(
								"<html><b>"+key+"</b></html>"
								,true);
						getTreeNodeRoot().add(n1);
						for(Object v2:o)
							{
							DefaultMutableTreeNode n2=new DefaultMutableTreeNode(
									String.valueOf(v2),
									false
									);
	
							n1.add(n2);
							}
						}
					}
				
				}
			}
		this.fireTreeStructureChanged();
		}
	private List<String> getCSQCols(VCFHeader header)
			{
			VCFInfoHeaderLine ihl=header.getInfoHeaderLine("CSQ");
			if(ihl==null) return Collections.emptyList();
			String description=ihl.getDescription();
			String chunck=" Format:";
			int i=description.indexOf(chunck);
			if(i==-1)
				{
				LOG.warn("Cannot find "+chunck+ " in "+description);
				return Collections.emptyList();
				}
			description=description.substring(i+chunck.length()).replaceAll("[ \'\\.\\(\\)]+","").trim();
			String tokens[]=pipeVep.split(description);
			ArrayList<String> L=new ArrayList<String>(tokens.length);
			for(String s:tokens)
				{
				if(s.trim().isEmpty()) continue;
				L.add(s);
				}
			return L;
			}
	private List<String> getEFFCols(VCFHeader header)
		{
		VCFInfoHeaderLine ihl=header.getInfoHeaderLine("EFF");
		if(ihl==null) return Collections.emptyList();
		String description=ihl.getDescription();
		String chunck="Format:";
		int i=description.indexOf(chunck);
		if(i==-1)
			{
			LOG.warn("Cannot find "+chunck+ " in "+description);
			return Collections.emptyList();
			}
		description=description.substring(i+chunck.length()).replace('(','|').replaceAll("[ \'\\.)\\[\\]]+","").trim();
		String tokens[]=pipeSnpEff.split(description);
		ArrayList<String> L=new ArrayList<String>(tokens.length);
		for(String s:tokens)
			{
			if(s.trim().isEmpty()) continue;
			L.add(s);
			}
		return L;
		}
	
	
	public void fireTreeStructureChanged()
		{
		fireTreeStructureChanged(this, new Object[]{getRoot()}, null, null);
		}	
	}
