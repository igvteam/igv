/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package util;

import org.junit.Ignore;
import org.junit.experimental.categories.Category;
import org.junit.internal.AssumptionViolatedException;
import org.junit.runner.Description;
import org.junit.runner.notification.RunNotifier;
import org.junit.runners.BlockJUnit4ClassRunner;
import org.junit.runners.model.FrameworkMethod;
import org.junit.runners.model.InitializationError;
import org.junit.runners.model.Statement;

import java.io.IOException;
import java.util.List;

/**
 * Custom test runner, for categorizing/ignoring tests based on input flags.
 *
 * Because we depend on network resources, we don't necessarily want
 * to fail tests if they aren't available. Can set at launch whether we just log and ignore IOExceptions
 *
 * We can also annotate test methods with @LongRunning and skip those if desired.
 * User: jacob
 * Date: 2013-Feb-27
 */
@Ignore
public class IGVTestRunner extends BlockJUnit4ClassRunner {

    public static final String IGNORE_IOEXCEPTIONS_PROPERTY = "ignore.ioexceptions";
    public static final String INCLUDE_LONGRUNNING_PROPERTY = "include.longrunning";

    static boolean ignoreIOExceptions = false;
    static boolean includeLongRunning = false;

    static{
        ignoreIOExceptions = Boolean.parseBoolean(System.getProperty(IGNORE_IOEXCEPTIONS_PROPERTY, "" + ignoreIOExceptions));
        includeLongRunning = Boolean.parseBoolean(System.getProperty(INCLUDE_LONGRUNNING_PROPERTY, "" + includeLongRunning));
    }

    public static boolean isLongRunning(FrameworkMethod member){
        boolean isLongRunning = false;
        Category categoriesAnn = member.getAnnotation(Category.class);
        if(categoriesAnn != null){
            Class[] categories = categoriesAnn.value();
            for(Class category: categories){
                if(category.equals(LongRunning.class)){
                    isLongRunning = true;
                    break;
                }
            }
        }
        return isLongRunning;
    }

    @Override
    protected void runChild(FrameworkMethod method, RunNotifier notifier) {
        boolean isLongRunning = isLongRunning(method);

        if(isLongRunning && !includeLongRunning){
            Description description = describeChild(method);
            notifier.fireTestIgnored(description);
        }else{
            super.runChild(method, notifier);
        }
    }

    @Override
    protected Statement methodInvoker(final FrameworkMethod method, Object test) {
        Statement statement =  super.methodInvoker(method, test);
        return new MethodInvoker(statement);
    }

    @Override
    protected void collectInitializationErrors(List<Throwable> errors) {
        super.collectInitializationErrors(errors);

        /**
         The purpose here is to prevent throwing errors
         if we include a class which isn't a test class, and thus has no test methods
         Usually we annotate these with @Ignore. But unless with annotate the class
         with RunWith it doesn't use this runner anyway. So it's pointless
         */
//        List<Throwable> removeErrors = new ArrayList<Throwable>(1);
//        for (Throwable error : errors) {
//            if (error.getMessage().contains("No runnable methods")) {
//                removeErrors.add(error);
//            }
//        }
//
//        for(Throwable removeError: removeErrors){
//            errors.remove(removeError);
//        }
    }

    @Override
    protected Object createTest() throws Exception {
        return super.createTest();
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
