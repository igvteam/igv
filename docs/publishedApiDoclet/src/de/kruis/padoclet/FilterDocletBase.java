/*
 *  PublishedApiDoclet - a filter proxy for any javadoc doclet
 *  
 *  Copyright (C) 2007,2010  Anselm Kruis <a.kruis@science-computing.de>
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

/*
 * This file is based on a public domain source found at
 * http://www.sixlegs.com/blog/java/exclude-javadoc-tag.html
 * 
 * Many thanks to Chris Nokleberg for this piece of code.
 */

package de.kruis.padoclet;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import com.sun.javadoc.DocErrorReporter;
import com.sun.javadoc.Doclet;
import com.sun.javadoc.LanguageVersion;
import com.sun.javadoc.RootDoc;

import de.kruis.padoclet.util.AbstractOption;
import de.kruis.padoclet.util.HalfDynamicProxy;
import de.kruis.padoclet.util.HalfDynamicProxy.MessageInterface;

/**
 * This class is a base class for javadoc filter doclets. 
 * 
 *  A filter doclet is plugged between the javadoc core and a second
 *  doclet, that is used to create the output. The filter doclet modifies 
 *  the information about the java code to be documented. It can hide 
 *  java items (as done by the {@link de.kruis.padoclet.PublishedApiDoclet})
 *  or change the information in any other way.
 *  
 *  <p>
 *  This class is not a doclet by itself. It is intended as a base class for 
 *  a doclet. The derived class has to implement the static methods required by 
 *  a doclet and to setup the options ({@link FilterDocletBase.Option#register(AbstractOption)})
 *  and the proxy table ({@link de.kruis.padoclet.util.HalfDynamicProxy#setProxyClassTable(Class[][])}).
 *  See {@link de.kruis.padoclet.PublishedApiDoclet} for an example.
 *  
 * @author kruis
 *
 */
public class FilterDocletBase implements MessageInterface {
    
    /**
     * The name of the system property, that contains the name of the
     * delegate doclet. If this system property is unset, the default 
     * doclet (<code>com.sun.tools.doclets.standard.Standard</code>) is used.
     */
    public static final String PAD_DELEGATE_DOCLET_SYSTEM_PROPERTY = "PublishedApiDoclet.delegate";
    
    
    /**
     * This class is used to perform a lazy instantiation of the
     * delegate / formating doclet class.
     * 
     * @author kruis
     * @pad.exclude 
     */
    private static class DH {
		/**
         * holds the delegate doclet
         */
        public static final Class<?> delegateDoclet;
        static {
            String classname = System.getProperty(PAD_DELEGATE_DOCLET_SYSTEM_PROPERTY);
            if (classname == null || classname.length() == 0)
            	classname = "com.sun.tools.doclets.standard.Standard";
            Class<?> clazz = null;
            try {
                clazz = Thread.currentThread().getContextClassLoader().loadClass(classname);
            } catch (Exception e) {
                e.printStackTrace();
            } 
            delegateDoclet = clazz;
        }
    }

    
    /**
     * holds the error reporter provided by the javadoc core
     */
    private DocErrorReporter errorReporter;
   
    /**
     * Create a new <code>FilterDocletBase</code> instance.
     * 
     * This constructor is <code>protected</code>, because this
     * constructor is intended to be called by subclasses only.
     */
    protected FilterDocletBase() {
    }
   

    /**
     * @return Returns the errorReporter provided by the doclet core.
     */
    public final DocErrorReporter getErrorReporter() {
        return errorReporter;
    }
    /**
     * @param errorReporter The errorReporter to set.
     */
    public final void setErrorReporter(DocErrorReporter errorReporter) {
        this.errorReporter = errorReporter;
    }


	/* (non-Javadoc)
	 * @see de.kruis.padoclet.HalfDynamicProxy.MessageInterface#recive(java.lang.String)
	 */
	public void emitMessage(String theMessage, int priority) {
		if(priority == MessageInterface.PRIORITY_DEBUG) {
			errorReporter.printNotice(theMessage);		
		} else if (priority == MessageInterface.PRIORITY_WARN) {
			errorReporter.printWarning(theMessage);
		} else {
			errorReporter.printError(theMessage);
		}
	}   
    
    /**
     * Invoke a static method on the delegate doclet.
     * 
     * This method performs a few security checks.
     * 
     * @param name name of the method to be invoked
     * @param par an array, that contains the method parameters
     * @return the return value of the invoked method.
     */
    private static Object delegateDocletInvoke(String name, Object[] par) {
        try{
        Method[] docletmethods = Doclet.class.getMethods();
        for(int i=0;i<docletmethods.length;i++) {
            if(! docletmethods[i].getName().equals(name)) 
                continue;
            int modifiers = docletmethods[i].getModifiers();
            if (! (Modifier.isStatic(modifiers) && Modifier.isPublic(modifiers)))
                continue;
            Class<?>[] partypes = docletmethods[i].getParameterTypes();
            if (partypes.length != (par!=null ? par.length : 0)) 
                continue;
            for(int j=0;j<partypes.length;j++) {
                if (! (par[j] == null || partypes[j].isInstance(par[j]))) 
                    continue;
            }
            // OK, we have the right method signature
            Method m = DH.delegateDoclet.getMethod(name, partypes);
            modifiers = m.getModifiers();
            if (! (Modifier.isStatic(modifiers) && Modifier.isPublic(modifiers)))
                throw new NoSuchMethodException("Method is not public static: "+m.toString());
            if (! docletmethods[i].getReturnType().isAssignableFrom(m.getReturnType()))
                throw new NoSuchMethodException("Method has incompatible return type: "+m.toString());
            try {
            	return m.invoke(null,par);
            }catch(InvocationTargetException e) {
            	Throwable targetException = e.getTargetException();
            	if (targetException instanceof Exception) {
            		throw (Exception) targetException;
            	} else {
            		throw new RuntimeException(targetException);
            	}
            }
        }
        throw new NoSuchMethodException(name);
        } catch(RuntimeException e) {
            throw e;
        } catch(Exception e) {
            throw new RuntimeException(e);
        }
    }
    
    
    /**
     * Option handling for doclets.
     * 
     * This class holds static methods and data about the options given 
     * to the doclet. Additionally instances oft the class represent
     * single options.
     * 
     * @author kruis
     *
     */
    protected static class Option extends AbstractOption{
    	
    	/**
    	 * all filter doclet options start with this string.
    	 */
    	public final static String namePrefix = "-pad";
    	
    	
    	/**
    	 * Create a new option, that has a value.
    	 * 
    	 * @param name the name
    	 * @param defaultValue the default value
    	 * @param isTag set to <code>true</code>, if the value of the option names
    	 * a tag.
    	 * @param description the description of the option
    	 */
    	public Option(String name, String defaultValue, boolean isTag, String description) {
    		super(name, namePrefix, defaultValue, isTag, description);
    	}
    	
    	/**
    	 * Create a new option, that has no value (boolean option).
    	 * 
    	 * @param name the name
    	 * @param description the description of the option.
    	 */
    	public Option(String name, String description) {
    		super(name,namePrefix,description);
    	}
    	
    	/**
    	 * holds a (sorted) map of all known options
    	 */
    	private static Map<String,AbstractOption> options = new TreeMap<String, AbstractOption>();

    	/**
    	 * Register an option.
    	 * 
    	 * This method is intended to be used from a static initializer. It let
    	 * you define the set of possible option.
    	 * 
    	 * @param option the option to register.
    	 */
    	public static void register(AbstractOption option) {
    		register(option,options);
    	}
    	
    	/**
    	 * Get an option by name.
    	 * 
    	 * @param name the name of the option.
    	 * @return the option object. Returns <code>null</code>, if no option with the given  
    	 * name was registered. 
    	 */
    	public static AbstractOption get(String name) {
    		return get(name, options);
    	}
    	/**
    	 * Get a string made from the descriptions of all registered options.
    	 * 
    	 * @return the compiled descriptions.
    	 */
    	public static String getDescriptions() {
    		return getDescriptions(options);
    	}
    	
    	/**
    	 * Get all tags
    	 * 
    	 * @return a set containing all tag names, that is the values of all
    	 * options where the property <code>isTag</code> is set.
    	 */
    	public static Set<String> getTags() {
    		return getTags(options);
    	}
    	
    	
    	/**
    	 * Get the number of parameters this option takes.
    	 * 
    	 * @param name the name of the option.
    	 * @return 1, if the option takes no parameters, 2, if the option takes a parameter. If the option is unknown, return 0.
    	 */
    	public static int optionLength(String name) {
    		return optionLength(name, options);
    	}
    	/**
    	 * Initialize the option values.
    	 * 
    	 * @param docletoptions the options as provided by the javadoc core.
    	 * @see Doclet#validOptions(java.lang.String[][], com.sun.javadoc.DocErrorReporter)
    	 * @see RootDoc#options()
    	 */
    	public static void initOptions(String [][]docletoptions) {
    		initOptions(docletoptions, options);
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
    	 * @throws Throwable
    	 */
    	public static void initJavaBeanProperties(Object bean) throws Throwable {
    		initJavaBeanProperties(bean, options);
    	}
    	
    	// register default options
    	static {
    		register(new Option("NoTagOptions","Do not add -tag-options to the option list of the formating doclet."+LI
    				+"Use this option, if your formating doclet doesn't understand the \"-tag <tagname>:X\" option."));
    		register(new Option("Help","Show this help message."));
    	}
    }
    
   
    
    /**
     * Filter the command line options seen by the formating Doclet.
     * 
     * Remove the options, that are processed by the filter doclet and 
     * add options to suppress warnings about unknown tags.
     * 
     * @param options the options as provided by the javadoc core
     * @return the filtered options
    	 * @see Doclet#validOptions(java.lang.String[][], com.sun.javadoc.DocErrorReporter)
    	 * @see RootDoc#options()
     */
    protected static String[][] filterOptions(String [][] options) {
        // filter our own options
        List<String[]> filteredOptions = new ArrayList<String[]>();
        for(int i=0; i<options.length; i++) {
        	if (Option.optionLength(options[i][0]) == 0)
        		filteredOptions.add(options[i]);
        }
        if((! Option.get("NoTagOptions").isSet()) && optionLengthHelper("-tag",FilterDocletBase.class.getName()) == 2) {
        	// the -tag option of the standard doclet seems to be supported
        	Iterator<String> iterator = Option.getTags().iterator();
        	while(iterator.hasNext()) {
        		filteredOptions.add(new String[] {"-tag", iterator.next()+":X"});
        	}
        }
        return (String[][]) filteredOptions.toArray(new String[filteredOptions.size()][]);
    }

    /**
     * Helper method to ease the implementation of the doclet <code>validOptions</code> method.
     * 
     * This method provides all you need in order to implement 
     * {@link Doclet#validOptions(java.lang.String[][], com.sun.javadoc.DocErrorReporter)}.
     * 
     * @param options the options
     * @param reporter the errorReporter
     * @param showHelp if <code>true</code>, the method will show an online help message
     * and return <code>false</code>.
     * @param className the name of the filter doclet
     * @return <code>true</code>, if all options are valid. Otherwise show a help message 
     * and return <code>false</code>.
     * @throws java.io.IOException
     * @see Doclet#validOptions(java.lang.String[][], com.sun.javadoc.DocErrorReporter)
     */
    protected static boolean validOptionsHelper(String[][] options,
            DocErrorReporter reporter, boolean showHelp, String className) throws java.io.IOException {
    	Option.initOptions(options);
    	AbstractOption helpOption = Option.get("Help");
    	if (helpOption != null && helpOption.isSet()) {
    		showHelp = true;
    	}
        if (!((Boolean) delegateDocletInvoke("validOptions", new Object[] {filterOptions(options),reporter})).booleanValue()) {
        	showHelp = true;
        }
        if (showHelp) {
        	reporter.printNotice(Option.LF+
        			FilterDocletBase.class.getName()+ " options:"+
        			Option.LF+Option.getDescriptions());
        	return false;
        }
        return true;
    }

    /**
     * Helper method used to implement the doclet <code>optionLength</code>
     * method.
     * 
     * @param option the name of the option
     * @param className the name of the filter doclet
     * @return the length of the option
     * @see Doclet#optionLength(java.lang.String)
     */
    protected static int optionLengthHelper(String option, String className) {
    	int length = Option.optionLength(option);
    	if (length > 0)
    		return length;
        length = ((Integer) delegateDocletInvoke("optionLength",new Object[] {option})).intValue();
        if ("-help".equals(option)) {
            System.out.println(Option.LF+"Provided by "+className+" doclet:"
                    +Option.LF+Option.getDescriptions());
        }
        return length;
    }

	/**
	 * Helper method used to implement the doclet <code>languageVersion</code>
	 * method.
	 * 
	 * See {@link PublishedApiDoclet#languageVersion()} for an example on how to
	 * use this method.
	 * 
	 * @return the language version supported by the delegate doclet.
	 */
    protected static LanguageVersion languageVersionHelper() {
    	return (LanguageVersion) delegateDocletInvoke("languageVersion", null);
    }
    
    
    /**
     * Helper method used to implement the doclet <code>start</code> method.
     * 
     * See {@link PublishedApiDoclet#start(RootDoc)} for an example on how to use this
     * method.
     * 
     * @param root the RootDoc object as provided by the javadoc core.
     * @param fd a newly created filter doclet object
     * @return see {@link Doclet#start(com.sun.javadoc.RootDoc)}.
     * @throws java.io.IOException
     */
    protected static boolean startHelper(RootDoc root, FilterDocletBase fd) throws java.io.IOException {
        // process our options
        Option.initOptions(root.options());
        try {
			Option.initJavaBeanProperties(fd);
		} catch (Throwable e) {
			e.printStackTrace();
			root.printError(e.toString());
			return false;
		}
        fd.setErrorReporter(root);
        RootDoc filteredRootDoc = (RootDoc) HalfDynamicProxy.getHDPProxy(root, RootDoc.class, HalfDynamicProxy.stateFactory(fd,fd));
        fd.preDelegateStartHook(filteredRootDoc);
        return ((Boolean) delegateDocletInvoke("start", new Object[]{ filteredRootDoc })).booleanValue();
    }
	
    /**
     * Hook method called prior to the start-method of the delegate doclet.
     * 
     * This implementation is an empty method, that does nothing.
     * Override this method, if you want to perform any operations within this hook.
     * 
	 * @param filteredRootDoc the filtered RootDoc.
	 */
	protected void preDelegateStartHook(RootDoc filteredRootDoc) {
	}

	/**
     * This class is the base of all the 
     * HalfDynamicProxy classes for the javadoc *Doc interfaces.
     * 
     * All <code>com.sun.javadoc.*Doc</code> interfaces extend the 
     * {@link Comparable} interface. Additionally instances of these interfaces
     * must be comparable by reference. Therefore this class contains a 
     * matching {@link #compareTo(Object)} implementation.
     * 
     * @author kruis
     *
     */
    public static class ComparableHandler extends HandlerBase  {
		/* (non-Javadoc)
		 * @see java.lang.Comparable#compareTo(java.lang.Object)
		 */
		@SuppressWarnings("unchecked")
		public int compareTo(Object o) {
			return ((Comparable<Object>)target).compareTo(unwrap(o));
		}		
	}
	
	/**
     * This class is the base of all the 
     * HalfDynamicProxy classes for the javadoc interfaces.
     * 
     * @author kruis
     *
     */
    public static class HandlerBase extends HalfDynamicProxy  {
		/**
		 * print a debug message.
		 * 
		 * @param message
		 */
		protected void debug(String message) {
			FilterDocletBase pad = (FilterDocletBase) getHDPStateUserObject();
			pad.errorReporter.printNotice(message);
		}
	}
	
	
	
}
