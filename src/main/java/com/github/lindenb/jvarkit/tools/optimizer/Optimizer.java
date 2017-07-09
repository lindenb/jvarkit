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
package com.github.lindenb.jvarkit.tools.optimizer;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.InMemoryCompiler;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.IOUtil;
import jdk.nashorn.internal.parser.JSONParser;

/**
BEGIN_DOC


END_DOC
*/

public class Optimizer extends Launcher
	{
	private static final Logger LOG = Logger.build(Optimizer.class).make();
	
	
	@Parameter(names={"-A","--all"},description="Run all possible combinations")
	private boolean run_all_combinations = false;

	
	private Random random = new Random(0L);
	private List<VariableParam> variableParams = new ArrayList<>();
	private enum Status {STATUS_CONTINUE,STATUS_STOP};
	private long generation = 0L;
	private Solution bestSolution = null;
	private Constructor<?> solutionConstructor = null;

	
	public static abstract class Solution implements Comparable<Solution>
		{
		protected long generation;
		protected final Map<String,Object> params;
		public Solution(final long generation,final Map<String,Object> params)
			{
			this.params = Collections.unmodifiableMap(params);
			this.generation = generation;
			}
		public abstract int execute();
		public void delete() {}
		@Override
		public String toString()
			{
			return "Generation "+this.generation+"\n"+
					this.params.keySet().stream().
						map(S-> " \""+S+"\"\t"+this.params.get(S)).collect(Collectors.joining("\n" ));
			}
		}
	
	abstract class VariableParam
		{
		private final String name;
		VariableParam(final String name) {
			this.name = name;
			}
		public String getName()
			{
			return name;
			}
		abstract Iterable<VariableValue> getAll();
		abstract VariableValue get();
		}
	
	class VariableValue
		{
		final VariableParam param;
		final Object value;
		VariableValue(VariableParam param,Object value) {
			this.param = param;
			this.value = value;
			}
		}
	
	class StringListParam extends VariableParam
		{
		final List<String> values;
		StringListParam(final String name,final Set<String> values)
			{
			super(name);
			this.values = new ArrayList<>(values);
			}
		
		@Override Iterable<VariableValue> getAll()
			{
			return values.stream().
					map(S->new VariableValue(StringListParam.this,S)).
					collect(Collectors.toList());
			}

		@Override VariableValue get()
			{
			final int idx = random.nextInt(this.values.size());
			return new VariableValue(StringListParam.this,values.get(idx));
			}
		}
	
	class IntRangeParam extends VariableParam
		implements Iterable<VariableValue>
		{
		int min=0;
		int max=1;
		int shift=1;
		IntRangeParam(final String name,int min,int max,int shift)
			{
			super(name);
			}
		@Override Iterable<VariableValue> getAll()
			{
			return this;
			}
		
		@Override VariableValue get()
			{
			int v= this.min + random.nextInt(this.max - this.min);
			return new VariableValue(this,v);
			}
		
		@Override
		public Iterator<VariableValue> iterator()
			{
			return new MyIterator();
			}
		private class MyIterator extends AbstractIterator<VariableValue>{
			int curr=IntRangeParam.this.min;
			@Override
			protected VariableValue advance()
				{
				final int v=this.curr;
				if(v>=IntRangeParam.this.max) return null;
				this.curr += IntRangeParam.this.shift;
				return new VariableValue(IntRangeParam.this,v);
				}
			}
		}
	
	class DoubleRangeParam extends VariableParam
		implements Iterable<VariableValue>

		{
		double min=0;
		double max=1;
		double shift=0.1;
		DoubleRangeParam(final String name,double min,double max,double shift) {
			super(name);
			}
		@Override Iterable<VariableValue> getAll()
			{
			return this;
			}
		
		@Override VariableValue get()
			{
			double v= this.min + random.nextDouble()*(max-min);
			return new VariableValue(this,v);
			}
		@Override
		public Iterator<VariableValue> iterator()
			{
			return new MyIterator();
			}
		private class MyIterator extends AbstractIterator<VariableValue>{
			double curr=DoubleRangeParam.this.min;
			@Override
			protected VariableValue advance()
				{
				final double v=this.curr;
				if(v>DoubleRangeParam.this.max) return null;
				this.curr += DoubleRangeParam.this.shift;
				return new VariableValue(DoubleRangeParam.this,v);
				}
			}

		}
	
	class BooleanParam extends VariableParam
		{
		BooleanParam(final String s) {
			super(s);
			}
		@Override Iterable<VariableValue> getAll()
			{
			return Arrays.asList(
					new VariableValue(this,Boolean.TRUE),
					new VariableValue(this,Boolean.FALSE)
				);
			}
		
		@Override VariableValue get()
			{
			return new VariableValue(this,random.nextBoolean());
			}
		}
	
	private Solution creatSolution(final Map<String,VariableValue> values)
		{
		try
			{
			return Solution.class.cast(this.solutionConstructor.newInstance(++this.generation,values));
			}
		catch (InstantiationException | IllegalAccessException | IllegalArgumentException | InvocationTargetException e)
			{
			throw new RuntimeException(e);
			}
		}
	
	private Status challenge(final Map<String,VariableValue> values)
		{
		this.generation++;
		LOG.info("Challenge["+generation+"]: " + values.values().stream().map(VV->"\""+VV.param.getName()+"\":"+VV.value).collect(Collectors.joining("; ")));
		
		Solution sol = creatSolution(values);
		try
			{
			int ret = sol.execute();
			if(ret!=0)
				{
				LOG.warn("solution "+sol+" failed");
				sol=null;
				}
			}
		catch(Throwable err)
			{
			LOG.warn("solution "+sol+" failed");
			LOG.warn(err);
			sol=null;
			}
		if(sol!=null)
			{
			if(this.bestSolution!=null &&  sol.compareTo(this.bestSolution)<0)
				{
				LOG.debug("Deleteing "+this.bestSolution);
				this.bestSolution.delete();
				}
			LOG.debug("Best is now "+ sol);
			this.bestSolution = sol;
			}
		return Status.STATUS_CONTINUE;
		}
	
	private Status runAllRecursive(final Map<String,VariableValue> values,int param_idx)
		{
		if(param_idx == this.variableParams.size())
			{
			return challenge(values);
			}
		else
			{
			final VariableParam param = this.variableParams.get(param_idx);
			for(final VariableValue variableValue : param.getAll()) {
				values.put(param.getName(),variableValue);
				Status status = runAllRecursive(values,param_idx+1);
				if(status.equals(Status.STATUS_STOP)) return status;
				}
			return Status.STATUS_CONTINUE;
			}
		}
	
	private void runAll()
		{
		final Map<String,VariableValue> values = new HashMap<>();
		runAllRecursive(values,0);
		}
	
	private void runRandom()
		{
		for(;;)
			{
			final Map<String,VariableValue> values = new HashMap<>(this.variableParams.size());
			for(final VariableParam param : this.variableParams) {
				values.put(param.getName(),param.get());
				}
			Status status =  challenge(values);
			if(status.equals(Status.STATUS_STOP)) break;
			}
		}
	
	
	private void readParams(final File jsonFile) throws IOException {
		IOUtil.assertFileIsWritable(jsonFile);
		FileReader r= new FileReader(jsonFile);
		JsonParser jsonParser=new JsonParser();
		final JsonElement root = jsonParser.parse(r);
		for(JsonElement pe:root.getAsJsonObject().get("params").getAsJsonArray())
			{
			JsonObject paramObj = pe.getAsJsonObject();
			final String name = paramObj.get("name").getAsString();
			if(this.variableParams.stream().filter(V->V.getName().equals(name)).findAny().isPresent())
				{
				throw new IOException("duplicate param name="+name);
				}
			final VariableParam variableParam;
			final String type = paramObj.get("type").getAsString();
			if(type.equals("strings"))
				{
				final Set<String> values = new LinkedHashSet<>();
				for(Iterator<JsonElement> it=pe.getAsJsonObject().get("values").getAsJsonArray().iterator();
						it.hasNext();
						)
					{
					values.add(it.next().getAsString());
					}
				variableParam = new StringListParam(name,values);
				}
			else if(type.equals("int"))
				{
				variableParam = new IntRangeParam(name,
						paramObj.get("min").getAsInt(),
						paramObj.get("max").getAsInt(),
						paramObj.get("shift").getAsInt()
						);
				}
			else if(type.equals("double"))
				{
				variableParam = new DoubleRangeParam(name,
						paramObj.get("min").getAsDouble(),
						paramObj.get("max").getAsDouble(),
						paramObj.get("shift").getAsDouble()
						);
				}
			else if(type.equals("bool") || type.equals("boolean"))
				{
				variableParam = new BooleanParam(name);
				}
			else
				{
				throw new IOException("undefine param type="+type);
				}
			this.variableParams.add(variableParam);
			}
		jsonReader.close();
		r.close();
	}
	
	@Override
	public int doWork(final List<String> args)
		{
		try
			{
			final String javaCode="";
			final String className="CustomClass" + this.random.nextInt();

			final StringWriter codestr= new StringWriter();
			final PrintWriter w = new PrintWriter(codestr);
			w.println("import java.util.*;");
			w.println("public static class "+className+" extends "+Solution.class.getName().replace("$", ".")+ "{");
			w.println(className+"(final long generation,final Map<String,Object> params) {");
			w.println("super(generation,params");
			w.println("}");
			w.println(javaCode);
			w.println("}");
			w.flush();
			
			
			
			
			
			final InMemoryCompiler inMemoryCompiler = new InMemoryCompiler();
			final Class<?> clazz = inMemoryCompiler.compileClass(className, javaCode);
			this.solutionConstructor =clazz.getConstructor(Long.class,Map.class);
			
			
			
			if(this.run_all_combinations)
				{
				this.runAll();
				}
			else
				{
				this.runRandom();
				}
			if(this.bestSolution==null) {
				LOG.error("no solution was found");
				return -1;
				}
			return 0;
			}
		catch (final Exception e)
			{
			LOG.error(e);
			return -1;
			}
		finally
			{
			}
		}
	
	public static void main(final String[] args)
		{
		new Optimizer().instanceMainWithExit(args);
		}

	}
