/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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

* 2015 moving to knime and adding predictions snpeff and VEP
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.vcffilterjs;


import java.io.File;
import java.io.FileReader;

import javax.script.Bindings;
import javax.script.Compilable;
import javax.script.CompiledScript;
import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;

import com.github.lindenb.jvarkit.util.cli.GetOpt;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter3;


/**
 * Author: Pierre Lindenbaum PhD. @yokofakun
 * Motivation : base class for javascript+vcf filters
 */
public abstract class AbstractVcfJavascript
	extends AbstractVCFFilter3
	{
	private CompiledScript  script=null;
	private ScriptEngine engine=null;

	/** expression in file */
	private File SCRIPT_FILE=null;
	/** expression in string */
	private String SCRIPT_EXPRESSION=null;

	protected AbstractVcfJavascript()
		{
		
		}
	
	protected CompiledScript getCompiledScript()
		{
		return this.script;
		}
	
	protected ScriptEngine getScriptEngine() {
		return this.engine;
		}
	
	public void setScriptExpression(String expression) {
		SCRIPT_EXPRESSION = expression;
		}
	
	public void setScriptFile(File src)
		{
		SCRIPT_FILE = src;
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -e (script) javascript expression.");
		out.println(" -f (script) javascript file.");
		out.println(" -o output file (default stdout).");
		super.printOptions(out);
		}
	
	@Override
	protected String getGetOptDefault() {
		return super.getGetOptDefault()+"e:f:o:";
		}
	
	@Override
	protected GetOptStatus handleOtherOptions(int c, GetOpt opt, String[] args)
		{
		switch(c)
			{
			case 'o':this.setOutputFile(new File(opt.getOptArg())); return GetOptStatus.OK;
			case 'e':this.setScriptExpression(opt.getOptArg());return GetOptStatus.OK;
			case 'f':this.setScriptFile(new File(opt.getOptArg()));return GetOptStatus.OK;
			default: return super.handleOtherOptions(c, opt, args);
			}
		}
	

	@Override
	public int initializeKnime()
		{
		if(SCRIPT_EXPRESSION==null && SCRIPT_FILE==null)
			{
			error("undefined script");
			return -1;
			}
		if(SCRIPT_EXPRESSION!=null && SCRIPT_FILE!=null)
			{
			error("both javascript file/expr are set");
			return -1;
			}
		ScriptEngineManager manager = new ScriptEngineManager();
		this.engine = manager.getEngineByName("js");
		if(this.engine==null)
			{
			error("The embedded 'javascript' engine is not available in java. Do you use the SUN/Oracle Java Runtime ?");
			return -1;
			}
		
		try
			{
			Compilable compilingEngine = (Compilable)this.engine;
			this.script = null;
			if(SCRIPT_FILE!=null)
				{
				info("Compiling "+SCRIPT_FILE);
				FileReader r=new FileReader(SCRIPT_FILE);
				this.script=compilingEngine.compile(r);
				r.close();
				}
			else
				{
				info("Compiling "+SCRIPT_EXPRESSION);
				this.script=compilingEngine.compile(SCRIPT_EXPRESSION);
				}
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		return 0;
		}


	@Override
	public void checkKnimeCancelled() {
		}


	protected boolean accept(Bindings bindings) throws ScriptException
		{
		Object result = getCompiledScript().eval(bindings);
		if(result==null) return false;;
		if(result instanceof Boolean)
			{
			if(Boolean.FALSE.equals(result)) return false;
			}
		else if(result instanceof Number)
			{
			if(((Number)result).intValue()!=1) return false;
			}
		else
			{
			warning("Script returned something that is not a boolean or a number:"+result.getClass());
			 return false;
			}
		return true;
		}
	
	
	@Override
	public void disposeKnime() {
		this.engine=null;
		this.SCRIPT_EXPRESSION=null;
		this.SCRIPT_FILE=null;
		super.disposeKnime();
		}

	}
