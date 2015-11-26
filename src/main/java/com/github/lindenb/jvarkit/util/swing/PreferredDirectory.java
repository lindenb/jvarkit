package com.github.lindenb.jvarkit.util.swing;

import java.io.File;
import java.util.prefs.Preferences;

public class PreferredDirectory {
	private static final String PREFKEY="jvarkit.util.swing.preferred.directory";

	private PreferredDirectory()
		{
		}
	public static File get()
		{
		return get(PreferredDirectory.class);
		}
	public static File get(Class<?> clazz)
		{
		return get(PREFKEY,clazz);
		}
	public static File get(final String key,final Class<?> clazz)
		{
		try {
			Preferences prefs  = Preferences.userNodeForPackage(clazz);
			final String dirS = prefs.get(key, null);
			if(dirS==null || dirS.trim().isEmpty()) return null;
			return new File(dirS);
		} catch (Exception e) {
			return null;
			}
		}
	public static void update(final File f)
		{
		update(PreferredDirectory.class,f);
		}
	
	public static void update(final Class<?> clazz,final File f)
		{
		update(PREFKEY,PreferredDirectory.class,f);
		}
	
	public static void update(final String key,final Class<?> clazz,final File f)
		{
		if(f!=null && f.isFile()) 
			{
			update(key,clazz,f.getParentFile());
			return;
			}
		if(f!=null && !f.isDirectory()) return;
		try {
			Preferences prefs  = Preferences.userNodeForPackage(clazz);
			if(prefs==null) return;
		
			if(f==null || f.getPath().trim().isEmpty())
				{
				prefs.remove(key);
				}
			else
				{
				prefs.put(key, f.getPath());
				}
			prefs.sync();
			prefs.flush();
			}
		catch(Exception err)
			{	
			}
		}
	}
