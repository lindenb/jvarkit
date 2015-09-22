package com.github.lindenb.jvarkit.util.command;

import java.util.List;

import javax.swing.Icon;

import org.apache.commons.cli.Option;

public interface CommandFactory extends CommandInfo
	{
	public Icon getIcon();
	public List<Option> getOptionList();
	public Command createCommand();
	public int instanceMain(String args[]);
	public void instanceMainWithExit(String args[]);
	}
