package com.github.lindenb.jvarkit.tools.springbatch;

import java.io.File;
import java.lang.reflect.Constructor;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.function.Function;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.springframework.batch.core.ExitStatus;
import org.springframework.batch.core.StepExecution;
import org.springframework.batch.core.annotation.BeforeStep;
import org.springframework.batch.item.ItemProcessor;

import com.github.lindenb.jvarkit.lang.InMemoryCompiler;
import com.github.lindenb.jvarkit.lang.JvarkitException;

import htsjdk.samtools.util.CloserUtil;

public class JdkFilterItemProcessor<I,O> implements ItemProcessor<I, O>{
	private static final Log LOG = LogFactory.getLog(JdkFilterItemProcessor.class);
	private static final Random rand = new Random(System.currentTimeMillis());
	private static final Map<String,Class<?>> CODE2CLASS = Collections.synchronizedMap(new HashMap<String,Class<?>>());
	private File scriptFile = null;
	private String scriptExpression = null;
	private Function transformer = null;
	
	public File getScriptFile() {
		return scriptFile;
		}
	
	public void setScriptFile(File scriptFile) {
		this.scriptFile = scriptFile;
		}
	
	public void setExpression(String scriptExpression) {
		this.scriptExpression = scriptExpression;
	}
	
	public String getExpression() {
		return scriptExpression;
		}
	
	@BeforeStep
	public void beforeStep(final StepExecution stepExecution) {
		try {
			InMemoryCompiler javac = new InMemoryCompiler();
			
			String javaCode = InMemoryCompiler.getTheSourceCode(getExpression(), getScriptFile());
			final Class<?> clazz ;
			
			if( CODE2CLASS.containsKey(javaCode))
				{
				if(LOG.isWarnEnabled()) LOG.warn("Code already compiled");
				clazz = CODE2CLASS.get(javaCode);
				}
			else
				{
				final String className = "FunctionImpl"+ rand.nextInt()+ System.currentTimeMillis();
				if(LOG.isInfoEnabled()) LOG.info("Compiling\n" + InMemoryCompiler.beautifyCode(javaCode));

				clazz = javac.compileClass(className,
					javaCode.replace("__CLASS__", className));
				CODE2CLASS.put(javaCode, clazz);
				}
			final Constructor<?> ctor = clazz.getConstructor();
			final Object instance = ctor.newInstance();
			if(!(instance instanceof Function)) {
				throw new JvarkitException.ScriptingError("Not an instance of Function");
				}
			this.transformer = (Function)instance;
			}
		catch(final Exception err)
			{
			throw new RuntimeException(err);
			}
		}
	
	public O process(final I arg0) throws Exception 
		{
		return (O)this.transformer.apply(arg0);
		}
	
    public ExitStatus afterStep(StepExecution stepExecution)
    {
    	CloserUtil.close(this.transformer);
    	this.transformer = null;
    	return stepExecution.getExitStatus();
    }
}
