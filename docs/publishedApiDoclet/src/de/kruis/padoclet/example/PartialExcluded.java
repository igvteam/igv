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

package de.kruis.padoclet.example;

/**
 * This class would be excluded from the javadoc documentation,
 * if there were no explicitly included members.
 * 
 * @author kruis
 * @pad.exclude sample
 */
public class PartialExcluded {

	/**
	 * Create a new instance
	 * 
	 * @pad.include
	 */
	public PartialExcluded() {
		super();
	}

	/**
	 * Do nasty things.
	 * 
	 * Note: this function might change in a future release.
	 */
	@Deprecated
	public void doNastyThings() {
		System.out.println("shit !!!");
	}
	
	
	/**
	 * Get a handle.
	 * 
	 * @param arg
	 * @pad.include sample
	 */
	@Annotation(string="this is a string", undocumented="another string")
	public static HandleObject getHandle(String arg) {
		return new HandleObject(arg);
	}

	/**
	 * Get the name of a {@link HandleObject}.
	 * 
	 * @param ho the handle object.
	 * @return the name
	 * @pad.include sample
	 */
	public String getHandleName(HandleObject ho) {
		return ho.getName();
	}
	
}
