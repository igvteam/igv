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


package de.kruis.padoclet.example.excludedPackage;

/**
 * This class is secret. 
 * 
 * The class does very nasty things. Therefore the class and
 * its package are excluded from the documentation.
 * 
 * @author kruis
 *
 */
public class Secret {

	/**
	 * Create a new instance.
	 */
	public Secret() {
		super();
	}

	/**
	 * Foo a little bit.
	 */
	public void foo() {
		System.out.println("foo");
	}
}
