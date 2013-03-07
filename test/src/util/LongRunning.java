/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package util;

import org.junit.Ignore;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

/**
 * Category for tests which run for an extremely long amount of time. The
 * definition of "long" is left purposefully ambiguous.
 *
 * TODO Right now we only support usage on methods, implementing on classes is more work but totally doable
 *
 * Example
 *
 * public class MyTest extends AbstractHeadlessTest{
 *
 *     @Category(LongRunning.class)
 *     @Test
 *     public void myLongTest() throws Exception{
 *          //Do stuff which takes a long time
 *     }
 *
 * }
 *
 *
 * User: jacob
 * Date: 2013-Mar-05
 */
@Ignore("Is an annotation, not a test")
@Target(ElementType.METHOD)
@Retention(RetentionPolicy.RUNTIME)
public @interface LongRunning {
}
