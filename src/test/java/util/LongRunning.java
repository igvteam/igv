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
