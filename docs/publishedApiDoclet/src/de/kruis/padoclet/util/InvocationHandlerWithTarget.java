/*
 *  PublishedApiDoclet - a filter proxy for any javadoc doclet
 *  
 *  Copyright (C) 2006, 2010  Anselm Kruis <a.kruis@science-computing.de>
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


/**
 * 
 */
package de.kruis.padoclet.util;

import java.lang.reflect.InvocationHandler;

/**
 * An invocationHandler with the two additional properties <i>target</i> 
 * and <i>state</i>. 
 * 
 * @author kruis
 *
 */
public interface InvocationHandlerWithTarget extends InvocationHandler {
	/**
	 * Initialize this invocation handler.
	 * 
	 * @param target object, for which this object is a proxy.
	 * @param state an arbitrary object used to hold state information
	 */
	void setupInvocationHandler(Object target, Object state);
	
	
	/**
	 * Get the target object, for which this invocation handler is a 
	 * proxy for.
	 * 
	 * @return the value of the target parameter of 
	 * {@link #setupInvocationHandler(Object, Object)}.
	 */
	Object getInvocationTarget();
}
