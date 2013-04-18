/*
 *  PublishedApiDoclet - a filter proxy for any javadoc doclet
 *  
 *  Copyright (C) 2007, 2010  Anselm Kruis <a.kruis@science-computing.de>
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */

package de.kruis.padoclet.util;

import java.beans.Introspector;
import java.beans.PropertyDescriptor;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;

/**
 * Option handling for doclets.
 * 
 * This class holds static methods and data about the options given 
 * to the doclet. Additionally instances of the class represent
 * single options.
 * 
 * @author kruis
 *
 */
public abstract class AbstractOption {
	/**
	 * Delimiters used to separate multiple alternative names for the same
	 * function.  
	 */
	public static final String TAG_DELIMITER = " \t\n\r\f,";
		
	/**
	 * Holds one line.separator character. 
	 */
	public final static String LF = System.getProperty("line.separator");
	/**
	 * Holds a line.separator and some spaces. Used to format option descriptions.
	 */
	public final static String LI = LF+"            ";
	/**
	 * the name of the option without the namePrefix.
	 */
	public final String name;
	/**
	 * the common name namePrefix of this option
	 */
	public final String namePrefix;
	/**
	 * the value of the option. For boolean options the values <code>"true"</code>
	 * and <code>"false"</code> are used.
	 */
	public String value;
	
	/**
	 * the default value for the option. This value is used, if the 
	 * option is not given on the command line. For boolean options the values <code>"true"</code>
	 * and <code>"false"</code> are used.
	 */
	public final String defaultValue;
	
	/**
	 * <code>true</code>, if the option has no value parameter. Otherwise, the 
	 * option takes one parameter.
	 */
	public final boolean isBoolean;
	/**
	 * <code>true</code>, if the value of the option names a javadoc tag. 
	 * This information is used to add <code>-tag tagname:X</code> options to 
	 * the command line of the formating doclet.
	 */
	public final boolean isTag;
	/**
	 * a description of the option. Used for the online help message.
	 */
	public final String description;
	
	/**
	 * Create a new option, that has a value.
	 * 
	 * @param name the name
	 * @param namePrefix the common namePrefix of the option name
	 * @param defaultValue the default value
	 * @param isTag set to <code>true</code>, if the value of the option names
	 * a tag.
	 * @param description the description of the option
	 */
	public AbstractOption(String name, String namePrefix, String defaultValue, boolean isTag, String description) {
		this.isBoolean = false;
		this.name = name;
		this.namePrefix = namePrefix;
		this.value = this.defaultValue = defaultValue;
		this.isTag = isTag;
		this.description = description;
	}
	
	/**
	 * Create a new option, that has no value (boolean option).
	 * 
	 * @param name the name
	 * @param namePrefix the common namePrefix of the option name
	 * @param description the description of the option.
	 */
	public AbstractOption(String name, String namePrefix, String description) {
		this.isBoolean = true;
		this.namePrefix = namePrefix;
		this.name = name;
		this.value = this.defaultValue = Boolean.FALSE.toString();
		this.isTag = false;
		this.description = description;
	}
	
	/**
	 * Is a boolean option given?
	 * 
	 * @return <code>true</code>, if a boolean option is set (to
	 * be exact, if the value of the option is <code>"true"</code>. Otherwise
	 * returns false.
	 */
	public boolean isSet() {
		return Boolean.valueOf(value).booleanValue();
	}
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
        String strTmp = namePrefix+name+(isBoolean?"":" <value>");  
		return strTmp 
			+("                                  ".substring(strTmp.length()))
			+"Default value is: '"+defaultValue+"'"+LI+description;
	}
	
	/**
	 * Register an option.
	 * 
	 * This method is intended to be used from a static initializer. It let
	 * you define the set of possible option.
	 * 
	 * @param option the option to register.
	 * @param options the options map
	 */
	protected static void register(AbstractOption option, Map<String,AbstractOption> options) {
		options.put(Introspector.decapitalize(option.name),option);
	}
	/**
	 * Get an option by name.
	 * 
	 * @param name the name of the option.
	 * @param options the options map
	 * @return the option object. Returns <code>null</code>, if no option with the given  
	 * name was registered. 
	 */
	protected static AbstractOption get(String name, Map<String, AbstractOption> options) {
		return options.get(Introspector.decapitalize(name));
	}
	/**
	 * Get a string made from the descriptions of all registered options.
	 * @param options the options map
	 * 
	 * @return the compiled descriptions.
	 */
	protected static String getDescriptions(Map<String, AbstractOption> options) {
		StringBuffer sb = new StringBuffer();
		Iterator<AbstractOption> iterator = options.values().iterator();
		while(iterator.hasNext()) {
			sb.append(iterator.next()).append(LF);
		}
		return sb.toString();
	}
	
	/**
	 * Get all tags
	 * @param options the options map
	 * 
	 * @return a set containing all tag names, that is the values of all
	 * options where the property <code>isTag</code> is set.
	 */
	protected static Set<String> getTags(Map<String, AbstractOption> options) {
		Set<String> set = new HashSet<String>();
		Iterator<AbstractOption> iterator = options.values().iterator();
		while(iterator.hasNext()) {
			AbstractOption o = iterator.next();
			if(o.isTag) {
				StringTokenizer tokenizer = new StringTokenizer(o.value,TAG_DELIMITER);
				while(tokenizer.hasMoreTokens()) {
					set.add(tokenizer.nextToken());
				}
			}
		}
		return set;
	}
	
	/**
	 * Get an option by its prefixed name.
	 * @param prefixedName the name of the option including the namePrefix.
	 * @param options the options map
	 * @return the option or <code>null</code>, if no matching option exists.
	 */
	private static AbstractOption getWithPrefix(String prefixedName, Map<String, AbstractOption> options) {
		if (prefixedName == null)
			return null;
		
		// This is a fairly slow implementation, but I don't consider option 
		// processing an very important issue.
		Iterator<Map.Entry<String, AbstractOption>> iterator = options.entrySet().iterator();
		while(iterator.hasNext()) {
			Map.Entry<String, AbstractOption> entry = iterator.next();
			String name = entry.getKey();
			AbstractOption o = entry.getValue();
			if (prefixedName.startsWith(o.namePrefix) 
					&& Introspector.decapitalize(prefixedName.substring(o.namePrefix.length())).equals(name))
				return o;
		}
		return null;
	}
	
	
	/**
	 * Get the number of parameters this option takes.
	 * 
	 * @param name the name of the option.
	 * @param options the options map
	 * @return 1, if the option takes no parameters, 2, if the option takes a parameter. If the option is unknown, return 0.
	 */
	protected static int optionLength(String name, Map<String, AbstractOption> options) {
		AbstractOption option = getWithPrefix(name, options);
		if (option == null)
			return 0;
		return option.isBoolean ? 1 : 2;
	}
	/**
	 * Initialize the option values.
	 * 
	 * @param options the options as provided by the javadoc core.
	 * @param optionsMap the options map
	 * @see com.sun.javadoc.Doclet#validOptions(java.lang.String[][], com.sun.javadoc.DocErrorReporter)
	 * @see com.sun.javadoc.RootDoc#options()
	 */
	protected static void initOptions(String [][]options, Map<String, AbstractOption> optionsMap) {
		for(int i=0;i<options.length;i++) {
	   		AbstractOption option = getWithPrefix(options[i][0], optionsMap);
	   	    if (option == null)
	   	    	continue;
	   	    if (option.isBoolean) {
	   	    	option.value = Boolean.toString(true);
	   	    } else {
	   	    	option.value = options[i][1];
	   	    }
		}
	}
	
	/**
	 * Assign option values to matching bean properties.
	 * 
	 * For each setable property of the Java bean, this method looks for an 
	 * option with the same name. If such an option exists, the property is set to the 
	 * value of the option. Currently only beans of the types <code>String</code>, <code>boolen</code> and 
	 * <code>int</code> are supported.
	 * 
	 * @param bean a java bean
	 * @param options the options map
	 * @throws Throwable
	 */
	protected static void initJavaBeanProperties(Object bean, Map<String, AbstractOption> options) throws Throwable {
		PropertyDescriptor[] pd = 
			Introspector.getBeanInfo(bean.getClass()).getPropertyDescriptors();
		for(int i=0;i<pd.length;i++) {
			Method wm = pd[i].getWriteMethod();
			if (wm == null)
				continue;
			String name = pd[i].getName();
			//System.err.println("Going to set Property: "+name);
			AbstractOption option = AbstractOption.get(name, options);
			if (option == null)
				continue;
			Class<?> propertyType = pd[i].getPropertyType();
			Object value = null;
			if (propertyType.isAssignableFrom(String.class)) {
				value = option.value;
			} else if (propertyType.isAssignableFrom(Integer.TYPE)) {
				value = new Integer(option.value);
			} else if (propertyType.isAssignableFrom(Boolean.TYPE)) {
				value = Boolean.valueOf(option.isSet());
			} else {
				continue;
			}
			try {
				wm.invoke(bean,new Object[]{ value});
    			//System.err.println("Done with Property: "+name+" new value is: "+value);
			}catch (InvocationTargetException e) {
				throw e.getTargetException();
			}
		}
	}
}