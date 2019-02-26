package com.github.lindenb.jvarkit.tools.tests;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

import org.testng.ITestNGListener;
import org.testng.TestListenerAdapter;
import org.testng.TestNG;

/*
http://testng.org/doc/documentation-main.html#running-testng-programmatically
*/
public class MiniTestNG {
	private static final Logger LOG= Logger.getLogger("minitestng");
	final ITestNGListener tla = new TestListenerAdapter();
	final TestNG testng = new TestNG();
	
	private int registerTestClass( final Set<Class<?>> classes,final Class<?> clazz) {
		if(classes.contains(clazz)) return 0;
		LOG.info("adding "+clazz);
		classes.add(clazz);
		final AlsoTest also = clazz.getAnnotation(AlsoTest.class);
		if(also!=null && also.value()!=null && also.value().length>0)
			{
			for(final Class<?> c2:also.value()) {
				if(registerTestClass(classes,c2)!=0) return -1;
				}
			}
		return 0;
		}
	
	private int registerTestClass( final Set<Class<?>> classes,final String className) {
		try {
			final Class<?> clazz = Class.forName(className);
			if(classes.contains(clazz)) return 0;
			return registerTestClass(classes,clazz);
			}
		catch (final Exception e) {
			e.printStackTrace();
			return -1;
			}
		}
	
	private int doWork(final List<String> args) {
		try {			
			String outputdir = null;
			int idx=0;
			while(idx<args.size()) {
				if(args.get(idx).equals("-d") && idx+1 < args.size() )
					{
					outputdir = args.get(idx+1);
					args.remove(idx);
					args.remove(idx);//twice
					}
				else
					{
					++idx;
					}
				}
			
			final Set<Class<?>> classes=new HashSet<>();
			args.stream().forEach(C->registerTestClass(classes, C));
			if(classes.isEmpty()) {
				LOG.warning("no classes defined");
				return -1;
			}
			
			
			testng.setTestClasses(classes.toArray(new Class<?>[classes.size()]));
			testng.addListener(tla);
			if(outputdir!=null) testng.setOutputDirectory(outputdir);
			testng.run(); 
			return 0;
		} catch(Throwable err) {
			err.printStackTrace();
			return -1;
		}
	}
	
	public static void main(final String[] args) {
		{
		System.exit(new MiniTestNG().doWork(new ArrayList<>(Arrays.asList(args))));
		}

	}

}
