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

import org.junit.internal.AssumptionViolatedException;
import org.junit.runners.BlockJUnit4ClassRunner;
import org.junit.runners.model.FrameworkMethod;
import org.junit.runners.model.InitializationError;
import org.junit.runners.model.Statement;

import java.io.IOException;

/**
 * Custom test runner. Because we depend on network resources, we don't necessarily want
 * to fail tests if they aren't available. Can set at launch whether we just log and ignore IOExceptions
 * User: jacob
 * Date: 2013-Feb-27
 */
public class IGVTestRunner extends BlockJUnit4ClassRunner{

    public static final String IGNORE_IOEXCEPTIONS_PROPERTY = "ignore.ioexceptions";
    static boolean ignoreIOExceptions = false;

    static{
        ignoreIOExceptions = Boolean.parseBoolean(System.getProperty(IGNORE_IOEXCEPTIONS_PROPERTY, "" + ignoreIOExceptions));
    }

    @Override
    protected Statement methodInvoker(FrameworkMethod method, Object test) {
        Statement statement =  super.methodInvoker(method, test);
        return new MethodInvoker(statement);
    }

    public IGVTestRunner(Class<?> klass) throws InitializationError {
        super(klass);
    }

    private static class MethodInvoker extends Statement{

        private Statement statement ;

        MethodInvoker(Statement statement){
            this.statement = statement;
        }

        @Override
        public void evaluate() throws Throwable {
            try{
                statement.evaluate();
            }catch (IOException e){
                if(ignoreIOExceptions){
                    e.printStackTrace();
                    throw new AssumptionViolatedException("IOException: " + e.getMessage());
                }else{
                    throw e;
                }
            }
        }
    }
}
