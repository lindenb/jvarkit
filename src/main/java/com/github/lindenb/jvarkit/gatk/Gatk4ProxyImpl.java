package com.github.lindenb.jvarkit.gatk;

import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineArgumentParser;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.*;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;


public class Gatk4ProxyImpl  extends org.broadinstitute.hellbender.Main implements Gatk4Proxy {
    private static final Logger LOG = LogManager.getLogger(Main.class);

    @Override
	public void execute(final List<String> argv) throws Exception {
		final String[] args = argv.toArray(new String[argv.size()]);
		LOG.info(getCommandLineName()+":Executing:  gatk "+ String.join(" ",argv));
		final CommandLineProgram program =
			this.setupConfigAndExtractProgram(args, 
				this.getPackageList(),
				this.getClassList(),
				this.getCommandLineName()
				);
	    final Object result = Main.runCommandLineProgram(program, args);
		if(result==null) return;
		if(Boolean.TRUE.equals(result)) return;
		LOG.warn("Returned "+ result.getClass());
		LOG.error("Result is "+ result);
		final Throwable err= (result instanceof Throwable?Throwable.class.cast(result):null);
		if(err!=null) {
			throw new RuntimeException(err);
			}
		else
			{
			throw new RuntimeException("Failure");
			}
		}
    
	@Override
	protected String getCommandLineName() {
		return this.getClass().getSimpleName();
		}

	}
