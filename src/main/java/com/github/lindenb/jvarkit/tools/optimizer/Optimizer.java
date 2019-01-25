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
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.OpenJdkCompiler;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Iso8601Date;

/**
BEGIN_DOC

This is a Genetic-Programming-like parameters optimizer.

use must provide the scafold of a java code that will override `Solution` and a JSON file describing the
parameters.

## An example of json config

```json

{
"params":[
	{
	"name":"param1",
	"type":"int",
	"min": 1,
	"max":10,
	"shift":1
	},
	{
	"name":"param2",
	"type":"double",
	"min": 0.01,
	"max": 0.1,
	"shift":0.01
	}
  ]
}
```


## The base class Solution

```java
public static abstract class Solution implements Comparable<Solution>
	{
	protected long generation =  -1L;
	protected final Map<String,Object> params;
	public Solution(final Map<String,Object> params)
		{
		this.params = Collections.unmodifiableMap(params);
		}
	// eval the result. Must be implemented by the user
	public abstract int execute() throws Exception;
	// delete any file associated to this solution 
	public void delete() {}
	//mate a Solution with another, returns null if mating is not possible
	public Solution mate(final Solution another) {
		return null;
	}
	
	@Override
	public boolean equals(final Object obj) {
		if(obj==this) return true;
		if(obj ==null || !(obj instanceof Solution)) return false;
		final Solution other=Solution.class.cast(obj);
		return this.params.equals(other.params);
		}
	
	@Override
	public int hashCode() {
		return this.params.hashCode();
		}
	}
```

## An example of custom class extending `Solution`

`__BASE__` will be replaced by the base class name `Solution`.
`__CLASS__` will be replaced by the current generated class name.

The user's code will be inserted in the following template:

```
 1  import java.util.*;
 2  import java.io.*;
 3  import java.util.stream.*;
 4  import java.util.function.*;
 5  import htsjdk.samtools.util.*;
 6  import htsjdk.variant.variantcontext.*;
 7  import htsjdk.variant.vcf.*;
 8  import javax.annotation.processing.Generated;
 9  @Generated(value="Optimizer",date="2017-07-10T11:20:07+0200")
10  public class __CLASS__ extends  __BASE__ {
11  public __CLASS__(final Map<String,Object> params) {
12  super(params);
13  }
14      // user's code starts here 
(...)  
93     // user's code ends here 
94  }

```

in __CLASS__ User must implement:

* 'compareTo' to compare two solutions
* 'execute' to compute the result with the current params. Returns '0' on success.
* 'delete' remove resources associated to this Solution.

END_DOC
*/
@Program(
		name="optimizer",
		description="Genetic-Programming-like parameters optimizer",
		keywords={"genetic-programming"}
		)
public class Optimizer extends Launcher
	{
	private static final Logger LOG = Logger.build(Optimizer.class).make();
	
	
	@Parameter(names={"-A","--all"},description="Run all possible combinations")
	private boolean run_all_combinations = false;
	@Parameter(names={"-c","--code","--source"},description="User's code.",required=true)
	private File useSourceCode = null;
	@Parameter(names={"-seed","--random"},description="Random seed. -1 == current time")
	private long  randomSeed = -1L;
	
	
	private Random random = new Random(0L);
	private List<VariableParam> variableParams = new ArrayList<>();
	private enum Status {STATUS_CONTINUE,STATUS_STOP};
	private long generation = 0L;
	private final List<Solution> bestSolutions = new ArrayList<>();
	private Constructor<?> solutionConstructor = null;
	

	
	public static abstract class Solution implements Comparable<Solution>
		{
		protected long generation=-1L;
		protected final Map<String,Object> params;
		public Solution(final Map<String,Object> params)
			{
			this.params = Collections.unmodifiableMap(params);
			}
		public abstract int execute() throws Exception;
		/** delete any file associated to this solution */
		public void delete() {}
		
		/** mate a Solution with another, returns null if mating is not possible*/
		public Solution mate(final Solution another) {
			return null;
		}
		
		@Override
		public boolean equals(final Object obj) {
			if(obj==this) return true;
			if(obj ==null || !(obj instanceof Solution)) return false;
			final Solution other=Solution.class.cast(obj);
			return this.params.equals(other.params);
			}
		
		@Override
		public int hashCode() {
			return this.params.hashCode();
			}
		
		@Override
		public String toString()
			{
			return "Generation "+this.generation+"\n"+
					this.params.keySet().stream().
						map(S-> " \""+S+"\" :"+this.params.get(S)).collect(Collectors.joining("\n" ));
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
			this.min=min;
			this.max=max;
			this.shift=shift;
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
			this.min=min;
			this.max=max;
			this.shift=shift;

			}
		@Override Iterable<VariableValue> getAll()
			{
			return this;
			}
		
		@Override VariableValue get()
			{
			double v= this.min + random.nextDouble()*(this.max-this.min);
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
			final Solution sol = Solution.class.cast(
					this.solutionConstructor.newInstance(values.entrySet().stream().
							collect(Collectors.toMap(E->E.getKey(),E->E.getValue().value))
						));
			
			sol.generation = ++this.generation;
			return sol;
			}
		catch (InstantiationException | IllegalAccessException | IllegalArgumentException | InvocationTargetException e)
			{
			throw new RuntimeException(e);
			}
		}
	
	private void addSolution(final Solution sol)
		{
		if(sol==null) return;
		if(this.bestSolutions.contains(sol)) return;
		this.bestSolutions.add(sol);
		Collections.sort(this.bestSolutions);
		if(this.bestSolutions.size()>10)
			{
			final Solution oldSol = this.bestSolutions.remove(this.bestSolutions.size()-1);
			LOG.debug("Deleting "+oldSol);
			oldSol.delete();
			}
		if(!this.bestSolutions.isEmpty())
			{
			LOG.debug("Best is\n"+ this.bestSolutions.get(0));
			}
		}
	
	private Status challenge(final Map<String,VariableValue> values)
		{
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
		catch(final Throwable err)
			{
			LOG.warn("solution "+sol+" failed");
			LOG.warn(err);
			sol.delete();
			sol=null;
			}
		if(sol!=null) addSolution(sol);
		
		// try mate
		if(!this.run_all_combinations && this.bestSolutions.size()>2)
			{
			int idx = 1 + this.random.nextInt(this.bestSolutions.size()-1);
			sol = this.bestSolutions.get(0).mate(this.bestSolutions.get(idx));
			if(sol!=null && !this.bestSolutions.contains(sol))
				{
				sol.generation = ++this.generation;
				try
					{
					int ret = sol.execute();
					if(ret!=0)
						{
						LOG.warn("solution "+sol+" failed");
						sol=null;
						}
					}
				catch(final Throwable err)
					{
					LOG.warn("solution "+sol+" failed");
					LOG.warn(err);
					sol.delete();
					sol=null;
					}
				if(sol!=null) addSolution(sol);
				}
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
	
	
	private int readParams(final File jsonFile) throws IOException {
		LOG.info("reading params in "+jsonFile);
		IOUtil.assertFileIsWritable(jsonFile);
		FileReader r= new FileReader(jsonFile);
		JsonParser jsonParser=new JsonParser();
		final JsonElement root = jsonParser.parse(r);
		if(!root.isJsonObject())
			{
			LOG.error("root is not an object");
			return -1;
			}
		if(!root.getAsJsonObject().has("params"))
			{
			LOG.error("root.params missing");
			return -1;
			}
		if(!root.getAsJsonObject().get("params").isJsonArray())
			{
			LOG.error("root.params is not an array");
			return -1;
			}
		if(root.getAsJsonObject().get("params").getAsJsonArray().size()==0)
			{
			LOG.error("root.params is  an empty array");
			return -1;
			}
		for(final JsonElement pe:root.getAsJsonObject().get("params").getAsJsonArray())
			{
			if(!pe.isJsonObject())
				{
				LOG.error("not an object:"+pe);
				return -1;
				}
			final JsonObject paramObj = pe.getAsJsonObject();
			if(root.getAsJsonObject().has("name"))
				{
				LOG.error("\"name\":missing in "+paramObj);
				return -1;
				}
			
			
			final String name = paramObj.get("name").getAsString();
			if(this.variableParams.stream().filter(V->V.getName().equals(name)).findAny().isPresent())
				{
				throw new IOException("duplicate param name="+name);
				}
			final VariableParam variableParam;
			if(root.getAsJsonObject().has("type"))
				{
				LOG.error("\"type\":missing in "+paramObj);
				return -1;
				}
			final String type = paramObj.get("type").getAsString();
			if(type.equals("strings"))
				{
				if(!paramObj.has("values"))
					{
					LOG.error("\"values\":missing in "+paramObj);
					return -1;
					}
				
				final Set<String> values = new LinkedHashSet<>();
				
				
				for(final Iterator<JsonElement> it= paramObj.get("values").getAsJsonArray().iterator();
						it.hasNext();
						)
					{
					values.add(it.next().getAsString());
					}
				if(values.isEmpty())
					{
					LOG.error("no values in "+paramObj);
					return -1;
					}
				variableParam = new StringListParam(name,values);
				}
			else if(type.equals("int") || type.equals("integer"))
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
		r.close();
		return 0;
		}
	
	@Override
	public int doWork(final List<String> args)
		{
		if(this.useSourceCode==null)
			{
			LOG.error("use source code is missing");
			return -1;
			}
		this.random =new Random(
				randomSeed==-1?
				System.currentTimeMillis():
				randomSeed
				);
		try
			{
			if(this.readParams( new File(super.oneAndOnlyOneFile(args)))!=0)
				{
				LOG.error("Cannot read params");
				return -1;
				}
			if(this.variableParams.isEmpty())
				{
				LOG.error("no params found");
				return -1;
				}
			
			final String className="CustomOptimizer" + Math.abs(this.random.nextInt());
			final StringWriter codestr= new StringWriter();
			final PrintWriter w = new PrintWriter(codestr);
			w.println("import java.util.*;");
			w.println("import java.io.*;");
			w.println("import java.util.stream.*;");
			w.println("import java.util.function.*;");
			w.println("import htsjdk.samtools.util.*;");
			w.println("import htsjdk.variant.variantcontext.*;");
			w.println("import htsjdk.variant.vcf.*;");
			w.println("import javax.annotation.processing.Generated;");
			w.println("@Generated(value=\""+Optimizer.class.getSimpleName()+"\",date=\""+ new Iso8601Date(new Date()) +"\")");
			w.println("public class "+className+" extends "+Solution.class.getName().replace("$", ".")+ "{");
			w.println("public "+className+"(final Map<String,Object> params) {");
			w.println("super(params);");
			w.println("}");
			w.println("   /** user's code starts here */");
			w.println(IOUtil.slurp(this.useSourceCode));
			w.println("   /** user's code ends here */");

			w.println("}");
			w.flush();
			
			final String code = codestr.toString().
					replaceAll("__BASE__", Solution.class.getName().replace("$", ".")).
					replaceAll("__CLASS__", className)
					;
			
			LOG.debug(" Compiling :\n" + OpenJdkCompiler.beautifyCode(code));

			
			final OpenJdkCompiler inMemoryCompiler = OpenJdkCompiler.getInstance();
			final Class<?> clazz = inMemoryCompiler.compileClass(className,code);
			this.solutionConstructor =clazz.getConstructor(Map.class);
			
			if(this.run_all_combinations)
				{
				this.runAll();
				}
			else
				{
				this.runRandom();
				}
			if(this.bestSolutions.isEmpty()) {
				LOG.error("no solution was found");
				return -1;
				}
			for(final Solution sol: this.bestSolutions)
				{
				sol.delete();
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
