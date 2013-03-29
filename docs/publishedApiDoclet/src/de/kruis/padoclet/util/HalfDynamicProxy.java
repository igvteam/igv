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

import java.lang.reflect.Array;
import java.lang.reflect.InvocationHandler;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.lang.reflect.Proxy;
import java.util.Map;
import java.util.WeakHashMap;


/**
 * This class is a quite universaly usable base class for a dynamic proxy.
 * 
 * @author kruis
 *
 */
public class HalfDynamicProxy implements InvocationHandlerWithTarget {
	/**
	 * Used to emit a message.
	 * 
	 * @author kruis
	 */
	public static interface MessageInterface {
		static int PRIORITY_ERROR = 0;

		static int PRIORITY_WARN = 1;

		static int PRIORITY_DEBUG = 2;

		void emitMessage(String theMessage, int priority);
	}

	/**
	 * Holds the state of a HalfDynamicProxy.
	 * 
	 * @author kruis
	 * 
	 */
	private static class HDPState {
		/**
		 * a default implementation of the {@link HalfDynamicProxy.MessageInterface}.
		 */
		private static HalfDynamicProxy.MessageInterface defaultReciver = new MessageInterface() {
			/*
			 * (non-Javadoc)
			 * 
			 * @see de.kruis.padoclet.HalfDynamicProxy.MessageInterface#emitMessage(java.lang.String,
			 *      int)
			 */
			public void emitMessage(String theMessage, int priority) {
				System.err.println(theMessage);
			}
		};

		/**
		 * Holds a user provided object. Not used by the HalfDynamicProxy class.
		 */
		private Object userState;

		/**
		 * holds the instance cache. This cache is required in order to build
		 * just a one proxy for each proxy target.
		 */
		private WeakHashMap<Object, Object> proxyCache;

		/**
		 * where to send messages to.
		 */
		private MessageInterface reciver;

		/**
		 * Create a new state object
		 * 
		 * @param userState
		 *            an arbitrary object provided by the caller. May be
		 *            <code>null</code>.
		 * @param reciver
		 *            where to send messages to. May be <code>null</code>.
		 */
		public HDPState(Object userState, MessageInterface reciver) {
			this.userState = userState;
			this.reciver = reciver != null ? reciver : defaultReciver;
			this.proxyCache = new WeakHashMap<Object, Object>();
		}

		@SuppressWarnings("unused")
		public void debug(String message) {
			this.reciver.emitMessage(message, MessageInterface.PRIORITY_DEBUG);
		}

		@SuppressWarnings("unused")
		public void error(String message) {
			this.reciver.emitMessage(message, MessageInterface.PRIORITY_ERROR);
		}
	}

	/**
	 * Holds a mapping table, that determinates the invocation handler class for
	 * a proxy target class.
	 */
	private static Class<?>[][] proxyClassTable;

	/**
	 * holds the target.
	 */
	protected Object target;

	/**
	 * holds the state
	 */
	private HalfDynamicProxy.HDPState state;

	/**
	 * Create a new state object.
	 * 
	 * @param userState
	 *            an arbitrary object provided by the caller. May be
	 *            <code>null</code>.
	 * @param reciver
	 *            where to send messages to. May be <code>null</code>.
	 * @return a new initialized HDPState object
	 */
	public static HDPState stateFactory(Object userState,
			MessageInterface reciver) {
		return new HDPState(userState, reciver);
	}

	/**
	 * @return Returns the proxyClassTable.
	 */
	public static Class<?>[][] getProxyClassTable() {
		return proxyClassTable;
	}

	/**
	 * Set the proxy class table.
	 * 
	 * The proxyClassTable is used to select the
	 * {@link InvocationHandlerWithTarget}, to be used for a given proxy target
	 * class. Each entry of the table is a <code>Class[]</code>, that
	 * contains a class implementing {@link InvocationHandlerWithTarget} as its
	 * first entry and zero or more interfaces as additional entries, that must
	 * be implemented by the object to create a proxy for. The first matching
	 * invocation handler is used.
	 * 
	 * @param proxyClassTable
	 *            The proxyClassTable to set.
	 */
	public static void setProxyClassTable(Class<?>[][] proxyClassTable) {
		HalfDynamicProxy.proxyClassTable = proxyClassTable;
	}

	/* (non-Javadoc)
	 * @see de.kruis.padoclet.InvocationHandlerWithTarget#setupInvocationHandler(java.lang.Object, java.lang.Object)
	 */
	public final void setupInvocationHandler(Object target, Object state) {
		this.target = target;
		this.state = (HDPState) state;
	}

	/**
	 * @return Returns the invovation target.
	 */
	public final Object getInvocationTarget() {
		return target;
	}

	/**
	 * Get the dynamic proxy for this object.
	 * @return an instance of {@link Proxy}.
	 */
	public Proxy dynamicProxyInstance() {
		return (Proxy) state.proxyCache.get(this.target);
	}

	/**
	 * Get the user object.
	 * 
	 * @return the user object, that was provided to the call of {@link #stateFactory(Object, HalfDynamicProxy.MessageInterface)}.
	 */
	protected Object getHDPStateUserObject() {
		return state.userState;
	}

	/**
	 * Get the proxy for an object using the same state as this proxy uses.
	 * 
	 * If <i>obj</i> is an array, this method 
	 * is called recursively on each element of the array.
	 * 
	 * @param obj the object to create a proxy for
	 * @param expect the expected type of the object,
	 * if <code>obj</code> is an array. If expect is {@link java.lang.Object} or 
	 * <code>null</code>, then the type of obj is used instead of expect.
	 * @return the proxy object, or obj itself.
	 * @see #getHDPProxy(Object, Class, HalfDynamicProxy.HDPState)
	 */
	protected Object getHDPProxy(Object obj, Class<?> expect) {
		return getHDPProxy(obj, expect, state);
	}

	/**
	 * Get the dynamic proxy for an object.
	 * 
	 *  If obj is an array, this method 
	 * is called recursively on each element of the array.
	 * No proxy is generated, if <code>obj</code> is <code>null</code> or primitive
	 * or, if <code>obj</code> is already a dynamic proxy, or, if
	 * no {@link InvocationHandlerWithTarget} implementation is found in the 
	 * {@link #proxyClassTable} (see {@link #setProxyClassTable(Class[][])}). 
	 *
	 * @param obj the object to create a proxy for.
	 * @param expect the expected type of the object,
	 * if <code>obj</code> is an array. If expect is {@link java.lang.Object} or 
	 * <code>null</code>, then the type of obj is used instead of expect.
	 * @param state the state object for the proxy
	 * @return the proxy object, or obj itself.
	 * @see #setProxyClassTable(Class[][])
	 */
	public static Object getHDPProxy(final Object obj, Class<?> expect, HDPState state) {
		if (obj == null) {
			return null;
		}

		// array handling
		if (obj instanceof Object[]) {
			Object[] arr = (Object[]) obj;
			if (null == expect || Object.class == expect) {
				// no explicit specification, use the type of the
				// array
				expect = obj.getClass();
			}
			Class<?> componentType = expect.getComponentType();
			boolean isProxy = true;
			for (int i = 0; i < arr.length; i++) {
				if (arr[i] != getHDPProxy(arr[i], componentType, state))
					isProxy = false;
			}
			if (isProxy)
				return obj;
			Object[] arr2 = (Object[]) Array.newInstance(componentType,
					arr.length);
			for (int i = 0; i < arr.length; i++) {
				arr2[i] = getHDPProxy(arr[i], componentType, state);
			}
			return arr2;
		}

		Class<?> cls = obj.getClass();
		if (Proxy.isProxyClass(cls)) {
			// is already a proxy
			return obj;
		}
		if (cls.isPrimitive() || obj instanceof String) {
			return obj;
		}
		// try to find an existing decorator
		Map<Object, Object> decoratorMap = state.proxyCache;
		synchronized (decoratorMap) {

			Object decorator = decoratorMap.get(obj);
			if (decorator != null) {
				return decorator;
			}
			// find the right class
			Class<?> invocationHandlerClass = null;
			for (int i = 0; i < proxyClassTable.length; i++) {
				// assert, that the row is valid
				if (proxyClassTable[i] == null
						|| proxyClassTable[i].length < 1
						|| !InvocationHandlerWithTarget.class
								.isAssignableFrom(proxyClassTable[i][0])) {
					throw new ClassCastException(
							"invalid proxy class at index " + i);
				}
				// loop over all required interfaces
				int j;
				for (j = 1; j < proxyClassTable[i].length; j++) {
					if (!proxyClassTable[i][j].isInstance(obj))
						break;
				}
				if (j >= proxyClassTable[i].length) {
					// we got it
					invocationHandlerClass = proxyClassTable[i][0];
					break;
				}
			}

			if (invocationHandlerClass == null) {
				// state.debug("no invocation Handler for object of class:
				// "+cls.getName());
				return obj;
			}
			InvocationHandlerWithTarget invokationHandler = null;
			try {
				invokationHandler = (InvocationHandlerWithTarget) invocationHandlerClass
						.newInstance();
				invokationHandler.setupInvocationHandler(obj, state);
			} catch (RuntimeException e) {
				throw e;
			} catch (Exception e) {
				e.printStackTrace();
				throw new RuntimeException(e);
			}

			Object proxy = Proxy.newProxyInstance(cls.getClassLoader(), cls
					.getInterfaces(), invokationHandler);
			// state.debug("created proxy: "+invocationHandlerClass.getName()+"
			// "+obj.toString());
			decoratorMap.put(obj, proxy);
			return proxy;
		}
	}

	/* (non-Javadoc)
	 * @see java.lang.reflect.InvocationHandler#invoke(java.lang.Object, java.lang.reflect.Method, java.lang.Object[])
	 */
	public Object invoke(Object proxy, Method method, Object[] args)
			throws Throwable {
		Method m = null;
		Object methodTarget = null;
		boolean isProxyRequired = false;
		try {
			m = this.getClass().getMethod(method.getName(),
					method.getParameterTypes());
			methodTarget = this;
			// state.debug("found replacement method "+m);
		} catch (NoSuchMethodException e) {
			m = method;
			methodTarget = this.target;
			isProxyRequired = true;
			// state.debug("no replacement method "+m);
		}
		Object result;
		try {
			result = m.invoke(methodTarget, args);
		} catch (InvocationTargetException e) {
			throw e.getTargetException();
		}
		if (isProxyRequired) {
			result = getHDPProxy(result, m.getReturnType());
		}
		return result;
	}

	/**
	 * Get the unwrapped proxy target.
	 * 
	 * @param proxy an object
	 * @return proxy, unless proxy is a dynamic proxy and its invocation handler
	 * implements the {@link InvocationHandlerWithTarget} interface. In 
	 * that case the target object will be returned.
	 */
	protected static Object unwrap(Object proxy) {
		if (proxy instanceof Proxy) {
			InvocationHandler invocationHandler = Proxy.getInvocationHandler(proxy);
			if (invocationHandler instanceof InvocationHandlerWithTarget) {
				 return ((InvocationHandlerWithTarget)invocationHandler)
				 .getInvocationTarget();
			}
		}
		return proxy;
	}

	/**
	 * This equals-method calls the equals 
	 * for the unwrapped proxy target objects.
	 * 
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	public boolean equals(Object obj) {
		return target.equals(unwrap(obj));
	}

	/**
	 * Calls the hashCode method of the proxy target. 
	 * 
	 * @see java.lang.Object#hashCode()
	 */
	public int hashCode() {
		return target.hashCode();
	}

	/** 
	 * Calls the toString method of the proxy target.
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
		return target.toString();
	}

}