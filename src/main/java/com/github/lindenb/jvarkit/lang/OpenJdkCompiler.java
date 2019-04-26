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

import java.io.File;
import java.io.InputStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;
import java.util.concurrent.TimeUnit;


import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;

public abstract class OpenJdkCompiler {
	private static final Logger LOG = Logger.build(OpenJdkCompiler.class).make();
	private static OpenJdkCompiler INSTANCE = null;
	private static final String JDK_PROPERTIES = "/META-INF/jdk.properties";
	public static OpenJdkCompiler getInstance() {
		if(INSTANCE==null) {
			synchronized (OpenJdkCompiler.class) {
				if(INSTANCE==null) {
					INSTANCE = new DefaultOpenJdkCompiler();
					}
				}
			}
		return INSTANCE;
		}
	
	
	public abstract Class<?> compileClass(final String className,final String javaCode);
	
	/** get full class name for the @Generated annotation. May be null. */
	public static String getGeneratedAnnotationClassName() {
		for(final String className: new String[]{
				"javax.annotation.Generated",
				"javax.annotation.processing.Generated"
				}) {
			try {
				final Class<?> c = Class.forName(className);
				if(c==null || !c.isAnnotation()) continue;
				return className;
				}
			catch(final ClassNotFoundException  err) {
				//ignore
				}
			}
		return null;
		}
	
	/** append line numbers to code */
	public static String beautifyCode(final String sourceCode)
		{
		final StringWriter codeWithLineNumber = new StringWriter();
		final String codeLines[] = sourceCode.split("[\n]");
		for(int nLine=0;nLine < codeLines.length;++nLine)
			{
			codeWithLineNumber.
				append(nLine==0?"":"\n").
				append(String.format("%10d  ",(nLine+1))+codeLines[nLine])
				;
			}
		return codeWithLineNumber.toString();
		}

	
	private static class DefaultOpenJdkCompiler
		extends OpenJdkCompiler
		{
		private Properties properties = null;
		DefaultOpenJdkCompiler() {
			}
		
		private Properties getProperties() {
			if(this.properties!=null) return properties;
			this.properties = new Properties();
			InputStream in = this.getClass().getResourceAsStream(JDK_PROPERTIES);
			if(in==null)
				{
				LOG.warn("cannot get resource: "+JDK_PROPERTIES);
				}
			else
				{
				try {
					this.properties.load(in);
					}
				catch(final Exception err) {
					throw new RuntimeIOException("Cannot read "+JDK_PROPERTIES,err);
					}
				finally
					{
					CloserUtil.close(in);
					}
				}
			return this.properties;
			}
		
		private String getJavacExe() {
			return getExecutable("javac");
			}
		private String getJarExe() {
			return getExecutable("jar");
			}
		private String getExecutable(final String name) {
			String s = getProperties().getProperty(name,null);
			if(!StringUtil.isBlank(s)) {
				final File exe = new File(s);
				if(exe.exists() && exe.isFile() && exe.canExecute()) return s;
				}
			
			try {
				String java_home_str = System.getenv("JAVA_HOME");
				if(!StringUtil.isBlank(java_home_str)) {
					File java_home = new File(java_home_str);
					if(java_home.exists() && java_home.isDirectory()) {
						File java_bin = new File(java_home,"bin");
						if(java_bin.exists() && java_bin.isDirectory()) {
							File exe = new File(java_bin,name);
							if(exe.exists() && exe.isFile() && exe.canExecute()) return exe.getPath();
							}
						}
					}
				else
					{
					LOG.warn("JAVA_HOME is not defined");
					}
				}
			catch(final SecurityException e) {
				// ignore
				}
			
			return name;
			}
		
		private void exec(final String definition,final List<String> cmd) {
			ProcessBuilder pb=new ProcessBuilder(cmd);
			try {
				Process proc = pb.
					redirectError(ProcessBuilder.Redirect.INHERIT).
					start();
				proc.waitFor(60L, TimeUnit.SECONDS);
				int ret = proc.exitValue();
				if(ret==0) return;
				}
			catch(Exception err)
				{
				LOG.error(err);
				}
			throw new RuntimeException("Cannot "+definition);
			}
		
		@Override
		public Class<?> compileClass(final String className,final String javaCode) {
			File javaSsrcDir = null;
			PrintWriter cw = null;
			File jarFile = null;
			try {
				//write source
				javaSsrcDir = IOUtil.createTempDir("jvarkit", ".tmp");
				IOUtil.assertDirectoryIsWritable(javaSsrcDir);
				final File javaFile = new File(javaSsrcDir,className+".java");
				cw = new PrintWriter(javaFile);
				cw.write(javaCode);
				cw.flush();
				cw.close();
				cw= null;
				//compile
				String classpath = getProperties().getProperty("classpath","");
                                boolean got_my_jar=false;
				//first try using this method before finding 'self', because jvarkit might be a symbolic link,
				// https://stackoverflow.com/questions/9399393/
				try
					{
					final java.security.CodeSource codeSource = getClass().getProtectionDomain().getCodeSource();
					final java.net.URL codeUrl = codeSource==null?null:codeSource.getLocation();
					final String myjar=codeUrl.toURI().getPath();
					if(myjar!=null && !myjar.isEmpty()) {
						classpath += File.pathSeparator+ myjar;
						got_my_jar=true;
						}
					}
				catch(final Throwable err)
					{
					LOG.debug("cannot get self jar using getProtectionDomain:"+err.getMessage());
					}

				if(!got_my_jar) {
					final String selfjar = getProperties().getProperty("self","");
					if(!selfjar.isEmpty()) {
						classpath += File.pathSeparator+ selfjar; 
						got_my_jar=true;
						}
					}

				if(!got_my_jar) {
					LOG.warn("cannot get self jar. Program might crash");
					}


				final List<String> cmd = new ArrayList<>();
				cmd.add(getJavacExe());
				cmd.add("-g");
				if(!classpath.isEmpty()) {
					cmd.add("-cp");
					cmd.add(classpath);
					}
				cmd.add("-d");
				cmd.add(javaSsrcDir.getPath());
				cmd.add("-sourcepath");
				cmd.add(javaSsrcDir.getPath());
				cmd.add(javaFile.getPath());
				exec("compile",cmd);
				
				//jar it
				jarFile = File.createTempFile("jvarkit", ".jar");
				jarFile.deleteOnExit();
				cmd.clear();
				cmd.add(getJarExe());
				cmd.add("cvf");
				cmd.add(jarFile.getPath());
				cmd.add("-C");
				cmd.add(javaSsrcDir.getPath());
				cmd.add(".");
				exec("jar",cmd);
				
				
				final URLClassLoader child = new URLClassLoader(
						new URL[] {jarFile.toURI().toURL()},
						this.getClass().getClassLoader()
						);
				final Class<?> compiledClass = Class.forName(className, true, child);
				return compiledClass;
				}
			catch(final Exception err) {
				throw new RuntimeException(err);
				}
			finally
				{
				CloserUtil.close(cw);
				if(javaSsrcDir!=null) IOUtil.deleteDirectoryTree(javaSsrcDir);
				}
			}
		}
	
	
}
