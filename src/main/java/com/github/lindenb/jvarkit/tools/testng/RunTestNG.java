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
package com.github.lindenb.jvarkit.tools.testng;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.testng.ITestNGListener;
import org.testng.TestListenerAdapter;

import org.testng.TestNG;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;

@Program(name="testng",
	description="run tests",
	keywords={"test","testng"},
	creationDate="202050512",
	modificationDate="202050512",
	generate_doc = false,
	jvarkit_hidden = true
	)
public class RunTestNG extends Launcher {
	private static final Logger LOG = Logger.of(RunTestNG.class);
	@Parameter(names={"-o","--output"},description= "output directory")
	private String outputDirectory = TestNG.DEFAULT_OUTPUTDIR;
	@Parameter(names={"--level"},description= "- the verbosity level (0 to 10 where 10 is most detailed). you can specify -1 and this will put TestNG in debug mode")
	private int verbose =5;
	@Parameter(names={"-i"},description= "skip/ignore missing classes")
	private boolean skip_missing_classes = false;

	private void loadClasses(BufferedReader br,final Set<Class<?>> m_test_classes) throws ClassNotFoundException,IOException {
		for(;;) {
			String line = br.readLine();
			if(line==null) break;
			if(line.startsWith("#") || StringUtils.isBlank(line)) continue;
			line=line.trim();
			try {
				final Class<?> c = Class.forName(line, false, this.getClass().getClassLoader());
				m_test_classes.add(c);
				}
			catch(ClassNotFoundException err) {
				if(skip_missing_classes) {
					LOG.warn("Class not found \""+line+"\" "+err.getMessage());
					continue;
					}
				throw err;
				}
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final Set<Class<?>> m_test_classes = new HashSet<>();
			
			final List<Path> paths  = IOUtils.unrollPaths(args);
			
			if(paths.isEmpty()) {
				try(BufferedReader r= new BufferedReader(new InputStreamReader(stdin()))) {
					loadClasses(r,m_test_classes);
					}
				}
			else
				{
				for(Path path:IOUtils.unrollPaths(args)) {
					try(BufferedReader r= Files.newBufferedReader(path)) {
						loadClasses(r,m_test_classes);
						}
					}
				}
			if(m_test_classes.isEmpty()) {
				LOG.warn("no class to test");
				return -1;
				}
			final ITestNGListener tla = new TestListenerAdapter();
			final TestNG testng = new TestNG();
			testng.setTestClasses(m_test_classes.stream().toArray(N->new Class<?>[N]));
			testng.addListener(tla);
			testng.setOutputDirectory(this.outputDirectory);
			testng.setVerbose(verbose);
			testng.run();
			return testng.hasFailure()?-1:0;
			}
		catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
    public static void main(String[] args) {
		new RunTestNG().instanceMainWithExit(args);
	}
    
}
