package com.github.lindenb.jvarkit.util.command;

import java.util.Collection;
import java.util.List;
import java.util.concurrent.Callable;

public interface Command extends CommandInfo,Callable<Collection<Throwable>>{
public CommandFactory getFactory();
public void setInput(List<String> nonOpts);
}
