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

/*
 * This file is based on a public domain source found at
 * http://www.sixlegs.com/blog/java/exclude-javadoc-tag.html
 * 
 * Many thanks to Chris Nokleberg for this piece of code.
 */

package de.kruis.padoclet;

import java.lang.reflect.Array;
import java.lang.reflect.Proxy;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.sun.javadoc.AnnotationDesc;
import com.sun.javadoc.AnnotationTypeDoc;
import com.sun.javadoc.AnnotationTypeElementDoc;
import com.sun.javadoc.AnnotationValue;
import com.sun.javadoc.ClassDoc;
import com.sun.javadoc.ConstructorDoc;
import com.sun.javadoc.Doc;
import com.sun.javadoc.DocErrorReporter;
import com.sun.javadoc.FieldDoc;
import com.sun.javadoc.LanguageVersion;
import com.sun.javadoc.MethodDoc;
import com.sun.javadoc.PackageDoc;
import com.sun.javadoc.Parameter;
import com.sun.javadoc.ProgramElementDoc;
import com.sun.javadoc.RootDoc;
import com.sun.javadoc.SourcePosition;
import com.sun.javadoc.Tag;
import com.sun.javadoc.Type;
import com.sun.tools.javadoc.Main;

import de.kruis.padoclet.util.AbstractOption;
import de.kruis.padoclet.util.HalfDynamicProxy;

/**
 * This class is a java 1.5 doclet, that is used as a filter between the javadoc
 * framework and another doclet, that produces some output.
 * 
 * The PublishedApiDoclet limits the packages, classes, fields, methods seen by
 * the second doclet based on exclude and include tags.
 * 
 * Technically this class is an application of the
 * {@link de.kruis.padoclet.FilterDocletBase} class.
 * 
 * @author kruis
 */
public class PublishedApiDoclet extends FilterDocletBase {
    
    /**
     * the name of the <i>exclude</i> tag. Default is <code>pad.exclude</code>.
     */
    private String excludeTag;
    /**
     * the exclude filter regular expression.
     */
    private Pattern excludeFilter;
    /**
     * the name of the <i>include</i> tag. Default is <code>pad.include</code>.
     */
    private String includeTag;
    /**
     * the include filter regular expression.
     */
    private Pattern includeFilter;
    /**
     * the name of the <i>forceInclude</i> tag. Default is <code>pad.forceInclude</code>.
     */
    private String forceIncludeTag;
    /**
     * the forceInclude filter regular expression.
     */
    private Pattern forceIncludeFilter;
    /**
     * the name of the <i>excludeChilds</i> tag. Default is <code>pad.excludeChilds</code>.
     */
    private String excludeChildsTag;
    /**
     * the excludeChilds filter regular expression.
     */
    private Pattern excludeChildsFilter;
    /**
     * the default mode. If <code>true</code>, everything is excluded by default. 
     * Default is <code>false</code> / include. 
     */
    private boolean defaultIsExclude;
    /**
     * the priority of the default setting. A positive value. Default is 1.
     */
    private int defaultPriority;
    /**
     * If <code>true</code>, retrieve all items from the javadoc core. Otherwise 
     * retrieve only the included subset.
     */
    private boolean disableJavadocFilter;
    /**
     * If <code>true</code>, don't call {@link Doc#isIncluded()}.
     */
    private boolean ignoreJavadocIsIncluded;
    
    /**
     * If <code>true</code>, do not filter enum constants
     */
    private boolean dontFilterEnumConstants;
    /**
     * If <code>true</code>, do not filter annotation elements
     */
    private boolean dontFilterAnnotationElements;
    /**
     * Create a new instance
     */
    private PublishedApiDoclet() {
    }
   
    /**
	 * @return Returns the defaultIsExclude.
	 */
	public final boolean isDefaultIsExclude() {
		return defaultIsExclude;
	}

	/**
	 * @param defaultIsExclude The defaultIsExclude to set.
	 */
	public final void setDefaultIsExclude(boolean defaultIsExclude) {
		this.defaultIsExclude = defaultIsExclude;
	}

	/**
	 * @return Returns the defaultPriority.
	 */
	public final int getDefaultPriority() {
		return defaultPriority;
	}

	/**
	 * @param defaultPriority The defaultPriority to set.
	 */
	public final void setDefaultPriority(int defaultPriority) {
		this.defaultPriority = defaultPriority;
	}

	/**
	 * @return Returns the disableJavadocFilter.
	 */
	public final boolean isDisableJavadocFilter() {
		return disableJavadocFilter;
	}

	/**
	 * @param disableJavadocFilter The disableJavadocFilter to set.
	 */
	public final void setDisableJavadocFilter(boolean disableJavadocFilter) {
		this.disableJavadocFilter = disableJavadocFilter;
	}

	/**
	 * @return Returns the excludeChildsFilter.
	 */
	public final String getExcludeChildsFilter() {
		return excludeChildsFilter.pattern();
	}
	public final Pattern getExcludeChildsFilterPat() {
		return excludeChildsFilter;
	}
	
	
	/**
	 * @param excludeChildsFilter The excludeChildsFilter to set.
	 */
	public final void setExcludeChildsFilter(String excludeChildsFilter) {
		this.excludeChildsFilter = Pattern.compile(excludeChildsFilter);
	}

	/**
	 * @return Returns the excludeChildsTag.
	 */
	public final String getExcludeChildsTag() {
		return excludeChildsTag;
	}

	/**
	 * @param excludeChildsTag The excludeChildsTag to set.
	 */
	public final void setExcludeChildsTag(String excludeChildsTag) {
		this.excludeChildsTag = excludeChildsTag;
	}

	/**
	 * @return Returns the excludeFilter.
	 */
	public final String getExcludeFilter() {
		return excludeFilter.pattern();
	}
	public final Pattern getExcludeFilterPat() {
		return excludeFilter;
	}
	
	/**
	 * @param excludeFilter The excludeFilter to set.
	 */
	public final void setExcludeFilter(String excludeFilter) {
		this.excludeFilter = Pattern.compile(excludeFilter);
	}

	/**
	 * @return Returns the forceIncludeFilter.
	 */
	public final String getForceIncludeFilter() {
		return forceIncludeFilter.pattern();
	}
	public final Pattern getForceIncludeFilterPat() {
		return forceIncludeFilter;
	}

	/**
	 * @param forceIncludeFilter The forceIncludeFilter to set.
	 */
	public final void setForceIncludeFilter(String forceIncludeFilter) {
		this.forceIncludeFilter = Pattern.compile(forceIncludeFilter);
	}

	/**
	 * @return Returns the forceIncludeTag.
	 */
	public final String getForceIncludeTag() {
		return forceIncludeTag;
	}

	/**
	 * @param forceIncludeTag The forceIncludeTag to set.
	 */
	public final void setForceIncludeTag(String forceIncludeTag) {
		this.forceIncludeTag = forceIncludeTag;
	}

	/**
	 * @return Returns the ignoreJavadocIsIncluded.
	 */
	public final boolean isIgnoreJavadocIsIncluded() {
		return ignoreJavadocIsIncluded;
	}

	/**
	 * @param ignoreJavadocIsIncluded The ignoreJavadocIsIncluded to set.
	 */
	public final void setIgnoreJavadocIsIncluded(boolean ignoreJavadocIsIncluded) {
		this.ignoreJavadocIsIncluded = ignoreJavadocIsIncluded;
	}

	/**
	 * @return Returns the includeFilter.
	 */
	public final String getIncludeFilter() {
		return includeFilter.pattern();
	}
	public final Pattern getIncludeFilterPat() {
		return includeFilter;
	}

	/**
	 * @param includeFilter The includeFilter to set.
	 */
	public final void setIncludeFilter(String includeFilter) {
		this.includeFilter = Pattern.compile(includeFilter);
	}

	/**
     * @return Returns the excludeTag.
     */
    public final String getExcludeTag() {
        return excludeTag;
    }
    /**
     * @param excludeTag The excludeTag to set.
     */
    public final void setExcludeTag(String excludeTag) {
        this.excludeTag = excludeTag;
    }
    /**
     * @return Returns the includeTag.
     */
    public final String getIncludeTag() {
        return includeTag;
    }
    /**
     * @param includeTag The includeTag to set.
     */
    public final void setIncludeTag(String includeTag) {
        this.includeTag = includeTag;
    }
	/**
	 * @return the dontFilterAnnotationElements
	 */
	public final boolean isDontFilterAnnotationElements() {
		return dontFilterAnnotationElements;
	}

	/**
	 * @param dontFilterAnnotationElements the dontFilterAnnotationElements to set
	 */
	public final void setDontFilterAnnotationElements(
			boolean dontFilterAnnotationElements) {
		this.dontFilterAnnotationElements = dontFilterAnnotationElements;
	}

	/**
	 * @return the dontFilterEnumConstants
	 */
	public final boolean isDontFilterEnumConstants() {
		return dontFilterEnumConstants;
	}

	/**
	 * @param dontFilterEnumConstants the dontFilterEnumConstants to set
	 */
	public final void setDontFilterEnumConstants(boolean dontFilterEnumConstants) {
		this.dontFilterEnumConstants = dontFilterEnumConstants;
	}

    
	// register the options. The option names must match the setable properties of 
    // the class
   	static {
    		Option.register(new Option("FilterDefault",".*",false,"Default regular expression for all filter options."+Option.LI
    				+"PublishedApiDoclet tags are ignored, unless the regular expression matches the text"+Option.LI
    				+"after the tag (that may be the empty string). The default regular expression \".*\" matches"+Option.LI
    				+"everything including the empty string."+Option.LI
    				+"  Every matched portion of the tag gets a priority number and the maximum of all these numbers"+Option.LI
    				+"is taken as the priority for this tag and language element. If a matched portion of the text"+Option.LI
    				+"contains a decimal number (e.g. \"foo123\") then the number is taken as the priority of the"+Option.LI
    				+"match. Otherwise the match priority is 1."+Option.LI
    				+"  The comparision of include- and exclude-priorities determinates, if the language element"+Option.LI
    				+"will be documented."));
    		Option.register(new Option("ExcludeTag","pad.exclude",true,"The name(s) of the exclude tag."+Option.LI
    				+"Use this option to redefine the tag. The empty string disables the tag."+Option.LI
    				+"Multiple names may be separated by white space or comma."));
    		Option.register(new Option("ExcludeFilter","<value of "+Option.namePrefix+"FilterDefault>",false,"The regular expression used to filter excludeTags."));
    		Option.register(new Option("IncludeTag","pad.include",true,"The name(s) of the include tag."+Option.LI
    				+"Use this option to redefine the tag. The empty string disables the tag."+Option.LI
    				+"Multiple names may be separated by white space or comma."));
    		Option.register(new Option("IncludeFilter","<value of "+Option.namePrefix+"FilterDefault>",false,"The regular expression used to filter includeTags."));
    		Option.register(new Option("ForceIncludeTag","pad.forceInclude",true,"The name(s) of the forceInclude tag. This tag is used,"+Option.LI
    				+"to include lenguage elements, that were excluded by the doclet itself (e.g. a private method)."+Option.LI
    				+"This tag has no priority."+Option.LI
    				+"  Use this option to redefine the tag. The empty string disables the tag."+Option.LI
    				+"Multiple names may be separated by white space or comma."));
    		Option.register(new Option("ForceIncludeFilter","<value of "+Option.namePrefix+"FilterDefault>",false,"The regular expression used to filter forceIncludeTags."));
    		Option.register(new Option("ExcludeChildsTag","pad.excludeChilds",true,"The name(s) of the excludeChilds tag. This tag is similar to the"+Option.LI
    				+"@pad.exclude tag, but it influences the exclude priority for childs of the taged element only."+Option.LI
    				+"The effective child exclude priority is the maximum of the @pad.exclude and the"+Option.LI
    				+"@pad.excludeChilds priorities."+Option.LI
    				+"  Use this option to redefine the tag. The empty string disables the tag."+Option.LI
    				+"Multiple names may be separated by white space or comma."));
    		Option.register(new Option("ExcludeChildsFilter","<value of "+Option.namePrefix+"FilterDefault>",false,"The regular expression used to filter excludeChildsTags."));
    		Option.register(new Option("DefaultIsExclude","Exclude all items except explicitly included ones. This option is usefull,"+Option.LI
    				+"if you want to include only a hand selected subset of your API."));
    		Option.register(new Option("DefaultPriority","1",false,"The priority of the default behaviour."));
    		Option.register(new Option("DisableJavadocFilter","Get unfiltered item collections from javadoc. You may need to"+Option.LI
    				+"use this option, if you use the @pad.forceInclude tag."));
    		Option.register(new Option("IgnoreJavadocIsIncluded","Do not call the javadoc isIncluded method."));
    		Option.register(new Option("DontFilterEnumConstants","Do not filter enum constants. (Enum constants are always included.)"));
    		Option.register(new Option("DontFilterAnnotationElements","Do not filter annotation elements. (Annotations elements are always included.)"));
		   	// make sure RefCheckDoclet is loaded and static initializers were execuded
		   	new RefCheckDoclet();
		   	AbstractOption option = RefCheckDoclet.Option.get(RefCheckDoclet.OPTION_WARN_ON);
		   	Option.register(new Option(option.name, "", option.isTag, option.description));
    	}
     
   	/**
   	 * Implements the doclet languageVersion method. 
   	 * 
   	 * This version of PA doclet only supports {@link LanguageVersion#JAVA_1_5}.
   	 * 
   	 * @return the supported language version
   	 */
   	public static LanguageVersion languageVersion() {
   		LanguageVersion lv = languageVersionHelper();
   		if (lv != LanguageVersion.JAVA_1_5) {
   			throw new RuntimeException("This version of "+PublishedApiDoclet.class.getName()+
   					" only supports language version "+LanguageVersion.JAVA_1_5+" but the "+
   					"delegate doclet uses language version "+lv);
   		}
   		return lv;
   	}

    /**
     * Implements the doclet validOptions method
     * 
     * @param options the options
     * @param reporter used to emit messages
     * @return <code>true</code>, is everything is ok, <code>false</code> otherwise.
     * @throws java.io.IOException
     * @see com.sun.javadoc.Doclet#validOptions(java.lang.String[][], com.sun.javadoc.DocErrorReporter)
     */
    public static boolean validOptions(String[][] options,
            DocErrorReporter reporter) throws java.io.IOException {
    	boolean showHelp = false;
    	
    	// init the options, because we need the DefaultPriority option
    	Option.initOptions(options);    	
    	try{
    		Integer.parseInt(Option.get("DefaultPriority").value);
    	} catch (NumberFormatException e) {
			reporter.printError("Option "+Option.namePrefix+"DefaultPriority"+" requires an integer argument" );
			showHelp = true;
    	}
    	// delegate most of the work to the helper method
        return validOptionsHelper(options,reporter,showHelp,PublishedApiDoclet.class.getName());
    }

    /**
     * The doclet optionLength method.
     * 
     * @param option the name of the option
     * @return the length or 0, if the option is unknown 
     * @see com.sun.javadoc.Doclet#optionLength(java.lang.String)
     */
    public static int optionLength(String option) {
    	// delegate the work to the helper method
        return optionLengthHelper(option,PublishedApiDoclet.class.getName());
    }

    /**
     * The doclet start method.
     * 
     * @param root the RootDoc object
     * @return <code>true</code>, if everything is ok, otherwise <code>false</code>.
     * @throws java.io.IOException
     * @see com.sun.javadoc.Doclet#start(com.sun.javadoc.RootDoc)
     */
    public static boolean start(RootDoc root) throws java.io.IOException {
        // process our options: set the default for the filter options
        Option.initOptions(root.options()); // first pass: we need FilterDefault
        Option.get("ExcludeFilter").value = 
        	Option.get("ExcludeChildsFilter").value =
        		Option.get("IncludeFilter").value =
        			Option.get("ForceIncludeFilter").value 
        			= Option.get("FilterDefault").value;
        
        // create the filter doclet instance
        FilterDocletBase fd = new PublishedApiDoclet();
        // delegate the work to the helper method
        return startHelper(root,fd);
    }

	/* (non-Javadoc)
	 * @see de.kruis.padoclet.FilterDocletBase#preDelegateStartHook(com.sun.javadoc.RootDoc)
	 */
	protected void preDelegateStartHook(RootDoc filteredRootDoc) {
		// the following lines are mor or less a copy of 
		// RefCheckDoclet#start()
		if (Option.get(RefCheckDoclet.OPTION_WARN_ON).value.length() > 0) {
			RefCheckDoclet rcd = new RefCheckDoclet();
			try {
				Option.initJavaBeanProperties(rcd);
			} catch (Throwable e) {
				e.printStackTrace();
				this.getErrorReporter().printError(e.toString());
				return;
			}
			rcd.setErrorReporter(this.getErrorReporter());
			rcd.check(filteredRootDoc);
		}
	}

    /**
     * A main method.
     * 
     * This method simply calls the doclet main method.
     * 
     * @param args the command line arguments
     * @see Main#execute(java.lang.String, java.lang.String, java.lang.String[])
     */
    public static void main(String[] args) {
        String name = PublishedApiDoclet.class.getName();
        Main.execute(name, name, args);
    }

    // initialize the proxy class table
	static { 
		HalfDynamicProxy.setProxyClassTable(
				new Class[][] {
				// 1st entry: handler class, 2nd and up: interfaces required to match
				// the interfaces are part of the com.sun.javadoc
				new Class[] { RootDocHandler.class,RootDoc.class},
				new Class[] { AnnotationTypeDocHandler.class, AnnotationTypeDoc.class},
				new Class[] { ClassDocHandler.class, ClassDoc.class},
				new Class[] { PackageDocHandler.class, PackageDoc.class},
				new Class[] { MethodDocHandler.class, MethodDoc.class},
				new Class[] { DocHandler.class, Doc.class},
				new Class[] { HandlerBase.class, AnnotationValue.class },
				new Class[] { AnnotationDescHandler.class, AnnotationDesc.class },
				new Class[] { HandlerBase.class, AnnotationDesc.ElementValuePair.class },
				new Class[] { ComparableHandler.class, Tag.class, Comparable.class},
				new Class[] { HandlerBase.class, Tag.class},
				new Class[] { HandlerBase.class, Type.class},
				new Class[] { HandlerBase.class, Parameter.class},
		}
		);
	}
	
	
	/**
	 * Proxy methods and state common to all {@link Doc} instances.
	 * 
	 * @author kruis
	 */
	public static class DocHandler extends ComparableHandler {


		/**
		 * The value of this field is only valid, if {@link #isIncludedValid} 
		 * is <code>true</code>.
		 * Then this field determinates, if the Doc instance is included.
		 */
		private boolean isIncluded;
		
		/**
		 * if <code>true</code>, the inclusion state is known and
		 * {@link #isIncluded} is valid
		 */
		private boolean isIncludedValid = false;
		
		/**
		 * Flag, that is used to avoid infinite recursion. May be obsolete.
		 */
		private boolean isCheckStarted = false;
		
		/**
		 * Cached inclusion priority for this node. 
		 * Only valid, if {@link #inclusionPriorityValid} is <code>true</code>.
		 */
		private int inclusionPriority = 0;
		/**
		 * Cached inclusion priority for childs of this node. 
		 * Only valid, if {@link #inclusionPriorityValid} is <code>true</code>.
		 */
		private int childInclusionPriority = 0;
		/**
		 * Indicates, if {@link #inclusionPriority} and {@link #childInclusionPriority} is
		 * valid.
		 */
		private boolean inclusionPriorityValid = false;
		
		/**
		 * Create a new instance.
		 */
		public DocHandler() {}
		
		/**
		 * Get the value of the boolean <i>included</i> property of the Doc 
		 * object, that is represented by this dynamic proxy instance.
		 * 
		 * This method is overridden by PublishedApiDoclet. It is the 
		 * switch used to turn documentation on or off.
		 * 
		 * @return return <code>true</code>, if this node (package, class, method, field) 
		 * shall be included in the documentation.
		 * 
		 * @see com.sun.javadoc.Doc#isIncluded()
		 */
		public boolean isIncluded() {
			// if the cached result is valid, there is no need to recompute the 
			//  result of this method.
			if (isIncludedValid)
				// the result is already known. Nothing to do.
				return isIncluded;
			
			// we probably do not need this check, but during development it proved 
			// to be useful. 
			if (isCheckStarted)
				throw new IllegalStateException("unexpected recursion detected");
			
			// start of try{ }finally
			try {
				// anti recursion flag.
				this.isCheckStarted = true;
				
				// get the state object and the target of this proxy
				PublishedApiDoclet pad = (PublishedApiDoclet) getHDPStateUserObject();
				Doc doc = (Doc)this.getInvocationTarget();
				
				// determinate, if the target is to be included in the 
				// documentation: Include it, if either 
				// 
				this.isIncluded = pad.isIgnoreJavadocIsIncluded() // the global flag says so
					|| doc.isIncluded()  // the javadoc framework says so
					|| // the "forceInclude"-tag of the target says so 
					   tagPriority(doc,pad.getForceIncludeTag(),pad.getForceIncludeFilterPat()) > 0;

			    // if we already know, that the target is not to be included,
			    // we are ready.
				if (! this.isIncluded)
					return false;

				// check the current target node and its parent nodes. 
				this.isIncluded = this.calcInclusionPriority(false) >= 0 ;
				
				// if the default is "exclude", check, if the
				// inclusion of just this target node is required,
				// because one of the child elements is included
				if (! this.isIncluded) {
					boolean inclusionRequired = false;
					if (doc instanceof PackageDoc) {
						// nested nodes are classes, therefore check, if the 
						// package contains any classes to be documented
						PackageDoc pd = (PackageDoc) dynamicProxyInstance();
						inclusionRequired = pd.allClasses().length > 0;
					} else if (doc instanceof ClassDoc) {
						// nested nodes are constructors, methods or fields.
						ClassDoc cd = (ClassDoc) dynamicProxyInstance();
						inclusionRequired = 
							cd.constructors().length > 0 ||
							cd.methods().length > 0 ||
							cd.fields().length > 0 ||
							cd.enumConstants().length > 0;
						if (! inclusionRequired && cd instanceof AnnotationTypeDoc) {
							AnnotationTypeDoc atd = (AnnotationTypeDoc) cd;
							inclusionRequired = atd.elements().length > 0;
						}
					}
					if (inclusionRequired) {
						debug("detected required inclusion of: "+doc.toString());
						this.isIncluded = true;
					}
				}	        
				return this.isIncluded;
			} finally {
				// mark the value of this.isIncluded as valid.
				this.isIncludedValid = true;
			}			
		}

		/**
		 * A constant used by {@link #tagPriority(Doc, String, Pattern)}.
		 */
		private static final Pattern digits = Pattern.compile(".*?(\\d+)");
		/**
		 * Get the priority from a tag.
		 * 
		 * If the doc-node does not contain a matching tag, the result is 
		 * 0. Otherwise the result is at least 1. If one or more matched portions of the 
		 * tag end with a decimal number the maximum of these numbers is returned. 
		 * 
		 * @param doc the document node (package, class, method, field, ...)
		 * @param tag the name of the tag
		 * @param filter a filter pattern. Only tags, that match this pattern are
		 * looked at.
		 * @return the priority asiciated with the tag: a non negative integer. 0 is the 
		 * minimum priority.
		 */
		private static int tagPriority(Doc doc, String tag, Pattern filter) {
			if (tag == null || tag.length() == 0)
				return 0;
			StringTokenizer tokenizer = new StringTokenizer(tag,Option.TAG_DELIMITER);
			int prio = 0;
			Tag[] tags;
			while(tokenizer.hasMoreTokens()) {
				tags = doc.tags(tokenizer.nextToken());
				int p;
				for(int i=0;i<tags.length;i++) {
					String text = tags[i].text();
					if (text == null)
						text=""; // I'm paranoid
					Matcher matcher = filter.matcher(text);
					while(matcher.find()) {
						// the tag matches, now get its priority
						p = 1;
						// now try to get the priority from the matched 
						matcher = digits.matcher(matcher.group());
						if (matcher.lookingAt()) {
							// the tag has an encoded priority
							p = Integer.parseInt(matcher.group(1));
						}
						// we need the maximum value
						if (prio < p)
							prio = p;
					}
				}
			}
			return prio;
		}
		
		/**
		 * Check, if this item shall be included.
		 * 
		 * This method checks this item and its parents, but not the childs.
		 * 
		 * @param iscallFromChild if <code>true</code> return the inclusion priority
		 * for a child item, otherwise return the inclusion priority for this item.
		 * @return an include (return value is positive) or exclude 
		 * (return value is negative) priority. 
		 */
		private int calcInclusionPriority(boolean iscallFromChild) {
			if (this.inclusionPriorityValid)
				return iscallFromChild ? this.childInclusionPriority : this.inclusionPriority;
			try {
				Doc doc = (Doc) dynamicProxyInstance();
				PublishedApiDoclet pad = (PublishedApiDoclet) getHDPStateUserObject();
				int includePriority = tagPriority(doc,pad.getIncludeTag(),pad.getIncludeFilterPat());
				int excludePriority = tagPriority(doc,pad.getExcludeTag(),pad.getExcludeFilterPat());
				int excludeChildsPriority = tagPriority(doc,pad.getExcludeChildsTag(),pad.getExcludeChildsFilterPat());
				if (excludePriority > excludeChildsPriority) {
					// the exclude childs priority must be at least as high 
					// as the exclude priority, because we do not want childs
					// to be included, and the local node included
					excludeChildsPriority = excludePriority;
				}
				
				int parentPriority = (pad.isDefaultIsExclude()?-1:1)*
						pad.getDefaultPriority();
				// ---- ask the parent, if possible ----
					
				Doc doc2 = null; // holds the parent
				// get the parent
				if (doc instanceof ProgramElementDoc) {
					ProgramElementDoc member = (ProgramElementDoc) doc;
					doc2 = member.containingClass();
					if (doc2 == null) {
						doc2 = member.containingPackage();
					}
				}
				// is parent valid?
				if (doc2 instanceof Proxy && 
						Proxy.getInvocationHandler(doc2) instanceof DocHandler) {
					// doc2 is the parent
					parentPriority = ((DocHandler)Proxy.getInvocationHandler(doc2)).calcInclusionPriority(true);
				}
				
				// now do the calculation for the local inclusion prio
				this.inclusionPriority = includePriority;
				if (excludePriority > includePriority) {
					this.inclusionPriority = -excludePriority;
				}
				if (Math.abs(parentPriority) > Math.abs(this.inclusionPriority)) {
					this.inclusionPriority = parentPriority;
				}
				// now do the calculation for the child inclusion prio
				this.childInclusionPriority = includePriority;
				if (excludeChildsPriority > includePriority) {
					this.childInclusionPriority = -excludeChildsPriority;
				}
				if (Math.abs(parentPriority) > Math.abs(this.childInclusionPriority)) {
					this.childInclusionPriority = parentPriority;
				}
				
				return iscallFromChild ? 
						this.childInclusionPriority : this.inclusionPriority;
			} finally {
				this.inclusionPriorityValid = true;
			}
		}
				
		/**
		 * Filter the content of a com.sun.javadoc.Doc-array, in order to remove entries not to be documented.
		 * 
		 * This method removes the entries from an array of Doc-items, 
		 * that were not to be documented. 
		 * 
		 * @param array the array of doc items.
		 * @param expect the class of the concrete array type, i.e. <code>com.sun.javadoc.PackageDoc[].class</code>, 
		 * if array contains PackageDoc items.
		 * @param doFilter if <code>false</code>, simply return a HD-proxy for array.
		 * @return a Doc[], containing HD-Proxies for the possibly filtered entries 
		 * of the given array.
		 */
		protected Doc[] filterDocArray(Doc[] array, Class<? extends Doc[]> expect, boolean doFilter) {
			if (!doFilter) {
				return (Doc[]) getHDPProxy(array, expect);
			}
			Class<?> componentType = expect.getComponentType();
			List<Doc> list = new ArrayList<Doc>(array.length);
			for (int i = 0; i < array.length; i++) {
				Doc entry = (Doc) getHDPProxy(array[i], componentType);
				if (entry != null && ! entry.isIncluded()){
					//debug("Array Excluding: "+entry.getClass().getName()+ " " + entry);
					continue;
				}
				list.add(entry);
			}
			return (Doc[]) list.toArray((Object[]) Array.newInstance(componentType,list.size()));
		}

		/**
		 * Determinate, if a com.sun.javadoc.Doc array shall be filtered.
		 * 
		 * @param filter if <code>false</code>, do not filter.
		 * @return returns <code>false</code>, if either the filter parameter is false or if 
		 * the filter is globally disabled.
		 */
		protected boolean isFilter(boolean filter) {
			if (! filter)
				return false;
			PublishedApiDoclet pad = (PublishedApiDoclet) getHDPStateUserObject();
			return ! pad.isDisableJavadocFilter();
		}

		/* (non-Javadoc)
		 * @see de.kruis.padoclet.FilterDocletBase.HandlerBase#debug(java.lang.String)
		 */
		protected void debug(String message) {
			SourcePosition position = ((Doc)(this.getInvocationTarget())).position();
			FilterDocletBase pad = (FilterDocletBase) getHDPStateUserObject();
			pad.getErrorReporter().printNotice(position,message);
		}

	}
	
	/**
	 * Proxy methods for the {@link MethodDoc} instance.
	 * 
	 * @author kruis
	 */
	public static class MethodDocHandler extends DocHandler {
		/**
		 * Create a new instance.
		 */
		public MethodDocHandler() {}
		
		/*
		 * @see MethodDoc#overrides(MethodDoc)
		 */
		public boolean overrides(MethodDoc meth) {
			return ((MethodDoc)target).overrides((MethodDoc) unwrap(meth));
		}
	}

	/**
	 * Proxy methods for the {@link RootDoc} instance.
	 * 
	 * @author kruis
	 */
	public static class RootDocHandler extends DocHandler {
		/**
		 * Create a new instance.
		 */
		public RootDocHandler() {}
		
		/* (non-Javadoc)
		 * @see com.sun.javadoc.RootDoc#classes()
		 */
		public ClassDoc[] classes() {
			return (ClassDoc[]) filterDocArray(((RootDoc)target).classes(),ClassDoc[].class,true);
		}

		/* (non-Javadoc)
		 * @see com.sun.javadoc.RootDoc#options()
		 */
		public String[][] options() {
			return filterOptions(((RootDoc)target).options());
		}
		
		/* (non-Javadoc)
		 * @see com.sun.javadoc.RootDoc#specifiedClasses()
		 */
		public ClassDoc[] specifiedClasses() {
			PublishedApiDoclet ds = (PublishedApiDoclet) getHDPStateUserObject();
			return (ClassDoc[]) filterDocArray(((RootDoc)target).specifiedClasses() ,ClassDoc[].class, 
					! ds.isDisableJavadocFilter());
		}

		/* (non-Javadoc)
		 * @see com.sun.javadoc.RootDoc#specifiedClasses()
		 */
		public PackageDoc[] specifiedPackages() {
			PublishedApiDoclet ds = (PublishedApiDoclet) getHDPStateUserObject();
			return (PackageDoc[]) filterDocArray(((RootDoc)target).specifiedPackages() ,PackageDoc[].class,
					! ds.isDisableJavadocFilter());
		}

	}

	/**
	 * Proxy methods for {@link ClassDoc} instances.
	 * 
	 * @author kruis
	 */
	public static class ClassDocHandler extends DocHandler {
		/**
		 * Create a new instance.
		 */
		public ClassDocHandler() {}
		
		/* (non-Javadoc)
		 * @see com.sun.javadoc.ClassDoc#constructors()
		 */
		public ConstructorDoc[] constructors() {
			return constructors(true);
		}
		/* (non-Javadoc)
		 * @see com.sun.javadoc.ClassDoc#constructors(boolean)
		 */
		public ConstructorDoc[] constructors(boolean filter) {
			return (ConstructorDoc[]) filterDocArray(((ClassDoc)target).constructors(isFilter(filter)) ,ConstructorDoc[] .class, filter);
		}
		/* (non-Javadoc)
		 * @see com.sun.javadoc.ClassDoc#fields()
		 */
		public FieldDoc[] fields() {
			return fields(true);
		}
		/* (non-Javadoc)
		 * @see com.sun.javadoc.ClassDoc#fields(boolean)
		 */
		public FieldDoc[] fields(boolean filter) {
			return (FieldDoc[]) filterDocArray(((ClassDoc)target).fields(isFilter(filter)) ,FieldDoc[] .class, filter);
		}
		/* (non-Javadoc)
		 * @see com.sun.javadoc.ClassDoc#innerClasses()
		 */
		public ClassDoc[] innerClasses() {
			return innerClasses(true);
		}
		/* (non-Javadoc)
		 * @see com.sun.javadoc.ClassDoc#innerClasses(boolean)
		 */
		public ClassDoc[] innerClasses(boolean filter) {
			return (ClassDoc[]) filterDocArray(((ClassDoc)target).innerClasses(isFilter(filter)) ,ClassDoc[] .class, filter);
		}
		/* (non-Javadoc)
		 * @see com.sun.javadoc.ClassDoc#methods()
		 */
		public MethodDoc[] methods() {
			return methods(true);
		}
		/* (non-Javadoc)
		 * @see com.sun.javadoc.ClassDoc#methods(boolean)
		 */
		public MethodDoc[] methods(boolean filter) {
			return (MethodDoc[]) filterDocArray(((ClassDoc)target).methods(isFilter(filter)) ,MethodDoc[] .class, filter);
		}
		/* (non-Javadoc)
		 * @see com.sun.javadoc.ClassDoc#subclassOf(com.sun.javadoc.ClassDoc)
		 */
		public boolean subclassOf(ClassDoc arg0) {
			return ((ClassDoc)target).subclassOf((ClassDoc) unwrap(arg0));
		}
		
		/*
		 * @see ClassDoc#enumConstants()
		 */
		public FieldDoc[] enumConstants() {
			PublishedApiDoclet pad = (PublishedApiDoclet) getHDPStateUserObject();
			return (FieldDoc[]) filterDocArray(((ClassDoc)target).enumConstants() ,FieldDoc[] .class, ! pad.isDontFilterEnumConstants());
		}
	}	
	
	/**
	 * Proxy methods for {@link ClassDoc} instances.
	 * 
	 * @author kruis
	 */
	public static class AnnotationTypeDocHandler extends ClassDocHandler {
		/**
		 * Create a new instance.
		 */
		public AnnotationTypeDocHandler() {}
		
		/*
		 * @see AnnotationTypeDoc#elements()
		 */
		public AnnotationTypeElementDoc[] elements() {
			PublishedApiDoclet pad = (PublishedApiDoclet) getHDPStateUserObject();
			return (AnnotationTypeElementDoc[]) filterDocArray(
					((AnnotationTypeDoc) target).elements(), AnnotationTypeElementDoc[].class,
					!pad.isDontFilterAnnotationElements());
		}
	}
	
	/**
	 * Proxy methods for {@link AnnotationDesc} instances.
	 * 
	 * @author kruis
	 */
	public static class AnnotationDescHandler extends HandlerBase {
		/**
		 * Create a new instance.
		 */
		public AnnotationDescHandler() {}
		
		/*
		 * @see AnnotationDesc#elementValues()
		 */
		public AnnotationDesc.ElementValuePair[] elementValues() {
			PublishedApiDoclet pad = (PublishedApiDoclet) getHDPStateUserObject();
			AnnotationDesc.ElementValuePair[] array = 
				((AnnotationDesc)target).elementValues();
			
			if (pad.isDontFilterAnnotationElements()) {
				return (AnnotationDesc.ElementValuePair[]) getHDPProxy(array, AnnotationDesc.ElementValuePair[].class);
			}
			List<AnnotationDesc.ElementValuePair> list = new ArrayList<AnnotationDesc.ElementValuePair>(array.length);
			for (int i = 0; i < array.length; i++) {
				AnnotationDesc.ElementValuePair pair = (AnnotationDesc.ElementValuePair) getHDPProxy(array[i], null);
				Doc entry = pair.element();
				if (entry != null && ! entry.isIncluded()){
					//debug("Array Excluding: "+entry.getClass().getName()+ " " + entry);
					continue;
				}
				list.add(pair);
			}
			return list.toArray((AnnotationDesc.ElementValuePair[]) Array.newInstance(AnnotationDesc.ElementValuePair.class,list.size()));
		}
	}
		
	/**
	 * Proxy methods for {@link PackageDoc} instances.
	 * 
	 * @author kruis
	 */
	public static class PackageDocHandler extends DocHandler {
		/**
		 * Create a new instance.
		 */
		public PackageDocHandler() {}
		
		/* (non-Javadoc)
		 * @see com.sun.javadoc.PackageDoc#allClasses()
		 */
		public ClassDoc[] allClasses() {
			return allClasses(true);
		}
		/* (non-Javadoc)
		 * @see com.sun.javadoc.PackageDoc#allClasses(boolean)
		 */
		public ClassDoc[] allClasses(boolean filter) {
			return (ClassDoc[]) filterDocArray(((PackageDoc)target).allClasses(isFilter(filter)) ,ClassDoc[] .class, filter);
		}
		/* (non-Javadoc)
		 * @see com.sun.javadoc.PackageDoc#errors()
		 */
		public ClassDoc[] errors() {
			return (ClassDoc[]) filterDocArray(((PackageDoc)target).errors() ,ClassDoc[] .class, true);
		}
		/* (non-Javadoc)
		 * @see com.sun.javadoc.PackageDoc#exceptions()
		 */
		public ClassDoc[] exceptions() {
			return (ClassDoc[]) filterDocArray(((PackageDoc)target).exceptions() ,ClassDoc[] .class, true);
		}
		/* (non-Javadoc)
		 * @see com.sun.javadoc.PackageDoc#interfaces()
		 */
		public ClassDoc[] interfaces() {
			return (ClassDoc[]) filterDocArray(((PackageDoc)target).interfaces() ,ClassDoc[] .class, true);
		}
		/* (non-Javadoc)
		 * @see com.sun.javadoc.PackageDoc#ordinaryClasses()
		 */
		public ClassDoc[] ordinaryClasses() {
			return (ClassDoc[]) filterDocArray(((PackageDoc)target).ordinaryClasses() ,ClassDoc[] .class, true);
		}
						
		public ClassDoc[] enums() {
			return (ClassDoc[]) filterDocArray(((PackageDoc)target).enums(), ClassDoc[].class, true);
		}
		
		public AnnotationTypeDoc[] annotationTypes() {
			return (AnnotationTypeDoc[]) filterDocArray(((PackageDoc)target).annotationTypes(), AnnotationTypeDoc[].class, true);
		}
		
	}

}
