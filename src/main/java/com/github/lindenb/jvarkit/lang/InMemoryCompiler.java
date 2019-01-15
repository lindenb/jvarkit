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
package com.github.lindenb.jvarkit.lang;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.StringWriter;
import java.net.URI;
import java.net.URL;
import java.util.Arrays;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.jar.Attributes;
import java.util.jar.Manifest;
import java.util.stream.Collectors;

import javax.tools.DiagnosticListener;
import javax.tools.FileObject;
import javax.tools.ForwardingJavaFileManager;
import javax.tools.JavaCompiler;
import javax.tools.JavaFileManager;
import javax.tools.JavaFileObject;
import javax.tools.SimpleJavaFileObject;
import javax.tools.StandardJavaFileManager;
import javax.tools.ToolProvider;

import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StringUtil;
/** 
 * In Memory Java Compiler
 * seen https://blog.nobel-joergensen.com/2008/07/16/using-eclipse-compiler-to-create-dynamic-java-objects-2/ 
 *
 */
@Deprecated //openjdk doesn't support javax.tools.JavaCompiler
public class InMemoryCompiler {
	private static final Logger LOG = Logger.build(InMemoryCompiler.class).make();
	/** custom class loader */
	private static class SpecialClassLoader extends ClassLoader {   
	    private final Map<String,MemoryByteCode> class2code = new HashMap<>();
	 
	    @Override
	    protected Class<?> findClass(final String name) throws ClassNotFoundException {       
	        MemoryByteCode mbc = this.class2code.get(name);       
	        if (mbc==null){           
	            mbc = class2code.get(name.replace(".","/"));           
	            if (mbc==null){               
	                return super.findClass(name);           
	            }       
	        }       
	        return super.defineClass(name, mbc.getBytes(), 0, mbc.getBytes().length);   
	    	}
	 
	    void addClass(final String name,final MemoryByteCode mbc) {       
	    	this.class2code.put(name, mbc);   
	    	}
		}
	
	/** custom SimpleJavaFileObject storing code in memory */ 
	private static class MemoryByteCode extends SimpleJavaFileObject {   
	    private ByteArrayOutputStream baos=null;   
	    public MemoryByteCode(final String name) {       
	        super(URI.create("byte:///" + name + ".class"), Kind.CLASS);   
	    	}   
	    
	    @Override
	    public CharSequence getCharContent(boolean ignoreEncodingErrors) {       
	        throw new IllegalStateException("Cannot get getCharContent");   
	    }   
	    
	    @Override
	    public OutputStream openOutputStream() {       
	        this.baos = new ByteArrayOutputStream();       
	        return this.baos;   
	    }   
	    @Override
	    public InputStream openInputStream() {       
	        throw new IllegalStateException();   
	    }   
	    
	    private byte[] getBytes() {       
	        return this.baos.toByteArray();   
	    }
	}
	 
	
	
	private static class MemorySource extends SimpleJavaFileObject {   
	    private final String src;   
	    public MemorySource(final String name, final String src) {       
	        super(URI.create("file:///" + name + ".java"), Kind.SOURCE);       
	        this.src = src;   
	    }   
	    @Override
	    public CharSequence getCharContent(boolean ignoreEncodingErrors) {       
	        return src;   
	    }   
	    @Override
	    public OutputStream openOutputStream() {       
	        throw new IllegalStateException();   
	    }   
	    @Override
	    public InputStream openInputStream() {       
	        return new ByteArrayInputStream(src.getBytes());   
	    }
	}
	
	private static class SpecialJavaFileManager extends ForwardingJavaFileManager<JavaFileManager> {   
	    private final SpecialClassLoader xcl;   
	    public SpecialJavaFileManager(final StandardJavaFileManager sjfm,final SpecialClassLoader xcl) {       
	        super(sjfm);       
	        this.xcl = xcl;   
	    }   
	    @Override
	    public JavaFileObject getJavaFileForOutput(final Location location,final String name, JavaFileObject.Kind kind, FileObject sibling) throws IOException {       
	        final MemoryByteCode mbc = new MemoryByteCode(name);       
	        xcl.addClass(name, mbc);       
	        return mbc;   
	    }
	    @Override
	    public ClassLoader getClassLoader(final Location location) {       
	        return xcl;   
	    }
	}
	 
	/** compile a new class */
	public Class<?> compileClass(final String className,final String javaCode)
		{
		 try{           
			final JavaCompiler javac = ToolProvider.getSystemJavaCompiler();
			if(javac==null)
				{
				throw new RuntimeException("ToolProvider.getSystemJavaCompiler() failed. Do you use a correct version of java ? Please check the version and avoid openJDK, use the java from Oracle.");
				}
			final StandardJavaFileManager sjfm = javac.getStandardFileManager(null, null, null);
			final SpecialClassLoader cl = new SpecialClassLoader();
			final SpecialJavaFileManager fileManager = new SpecialJavaFileManager(sjfm, cl);
			
			// https://stackoverflow.com/questions/1563909
			final Set<String> classpathcomponents = new LinkedHashSet<>(); 
			
			
			final Enumeration<URL> resources = getClass().getClassLoader().getResources("META-INF/MANIFEST.MF");
			while (resources.hasMoreElements()) {
				InputStream is = null;
			    try {
			        is = resources.nextElement().openStream();
			        final Manifest manifest = new Manifest(is);
			        final  Attributes attr = manifest.getMainAttributes();
			        final String cp=attr.getValue("Class-Path");
			        if(!StringUtil.isBlank(cp))
			        	{
			        	classpathcomponents.addAll( Arrays.stream(cp.split("[: ]")).
			        			filter(S->!S.trim().isEmpty()).
			        			collect(Collectors.toSet()));
			        	}
			    	} 
			    catch (final IOException err) {
			    	err.printStackTrace();
			    	}
			    finally
			    	{
			    	CloserUtil.close(is);
			    	}
				}
			
			final String java_class_path = System.getProperty("java.class.path");
			if(!StringUtil.isBlank(java_class_path)) {
	        	classpathcomponents.addAll( Arrays.stream(java_class_path.split("[: ]")).
	        			filter(S->!S.trim().isEmpty()).
	        			collect(Collectors.toSet()));
				}

			
			final List<String> options;
			if(!classpathcomponents.isEmpty())
				{
				final String classpathstr = String.join(":", classpathcomponents);
			    options = Arrays.asList("-classpath",classpathstr);
				}
			else
				{
				options = Collections.emptyList();
				}
			final List<JavaFileObject> compilationUnits = Arrays.asList(new MemorySource(className, javaCode));
			final DiagnosticListener<? super JavaFileObject> dianosticListener = null;
			final Iterable<String> classes = null;
			final StringWriter err = new StringWriter();
			final JavaCompiler.CompilationTask compile = javac.getTask(
					err, (JavaFileManager) fileManager, dianosticListener,
					options, classes, compilationUnits);
			if (compile.call()) {
				return cl.findClass(className);
				}
			else
				{
				LOG.warning(err.toString());
				throw new RuntimeException("compile.call() failed for custom class "+className);
				}
		} catch (final Exception e) {
			LOG.error(e);
			throw new RuntimeException("Cannot compile custom class "+className,e);
		}
	}
	
	/** append line numbers to code */
	public static String beautifyCode(final String sourceCode)
		{
		return StringUtils.beautifyCode(sourceCode);
		}
	
	/** get a code to compile in either expression or a file.
	 * Will throw an exception if BOTH expression and file are Both null or both not null */
	public static String getTheSourceCode(final String scriptExpr,final File scriptFile) throws IOException {
		if(scriptExpr==null && scriptFile==null) {
			throw new JvarkitException.ScriptingError("Both scriptExpr and scriptFile are undefined.");
			}
		else if(scriptExpr!=null && scriptFile!=null) {
			throw new JvarkitException.ScriptingError("Both scriptExpr and scriptFile are defined.");
			}
		else if(scriptExpr!=null)
			{
			return scriptExpr;
			}
		else
			{
			return IOUtil.slurp(scriptFile);
			}
		}
}
