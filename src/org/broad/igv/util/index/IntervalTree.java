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

package org.broad.igv.util.index;


/** An implementation of an interval tree, following the explanation.
 * from CLR.
 */


import org.apache.log4j.Logger;

import java.io.DataOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class IntervalTree {

    private static Logger logger = Logger.getLogger(IntervalTree.class);
    boolean immutable = false;
    Node root;
    Node NIL = Node.NIL;

    public IntervalTree() {
        this.root = NIL;
    }

    public IntervalTree(boolean immutable) {
        this.immutable = immutable;
        this.root = NIL;
    }

    public void insert(Interval interval) {
        if (immutable) {
            throw new java.lang.UnsupportedOperationException("Tree is immutable.  Inserts not allowed");
        }
        Node node = new Node(interval);
        insert(node);
    }


    // Returns all matches as a list of Intervals
    public List<Interval> findOverlapping(int start, int end) {

        Interval searchInterval = new Interval(start, end, 0);

        if (root().isNull()) {
            return Collections.emptyList();
        }

        List<Interval> results = new ArrayList();
        searchAll(searchInterval, root(), results);
        return results;
    }

    public String toString() {
        return root().toString();
    }

    private List<Interval> searchAll(Interval interval, Node node, List<Interval> results) {

        if (node.interval.overlaps(interval)) {
            results.add(node.interval);
        }

        if (!node.left.isNull() && node.left.max >= interval.getLow()) {
            searchAll(interval, node.left, results);
        }

        if (!node.right.isNull() && node.right.min <= interval.getHigh()) {
            searchAll(interval, node.right, results);
        }

        return results;
    }

    /**
     * Return all intervals in tree.
     * TODO: an iterator would be more effecient.
     * @return
     */
    public List<Interval> getIntervals() {
        if (root().isNull()) {
            return Collections.emptyList();
        }
        List<Interval> results = new ArrayList(size());
        getAll(root, results);
        return results;
    }

    private List<Interval> getAll(Node node, List<Interval> results) {

        results.add(node.interval);
        if (!node.left.isNull()) {
            getAll(node.left, results);
        }
        if (!node.right.isNull()) {
            getAll(node.right, results);
        }
        return results;
    }


    /**
     * Used for testing only.
     *
     * @param node
     * @return
     */
    private int getRealMax(Node node) {
        if (node.isNull())
            return Integer.MIN_VALUE;
        int leftMax = getRealMax(node.left);
        int rightMax = getRealMax(node.right);
        int nodeHigh = (node.interval).getHigh();

        int max1 = (leftMax > rightMax ? leftMax : rightMax);
        return (max1 > nodeHigh ? max1 : nodeHigh);
    }

    /**
     * Used for testing only
     *
     * @param node
     * @return
     */
    private int getRealMin(Node node) {
        if (node.isNull())
            return Integer.MAX_VALUE;

        int leftMin = getRealMin(node.left);
        int rightMin = getRealMin(node.right);
        int nodeLow = (node.interval).getLow();

        int min1 = (leftMin < rightMin ? leftMin : rightMin);
        return (min1 < nodeLow ? min1 : nodeLow);
    }


    private void insert(Node x) {
        assert (x != null);
        assert (!x.isNull());

        treeInsert(x);
        x.color = Node.RED;
        while (x != this.root && x.parent.color == Node.RED) {
            if (x.parent == x.parent.parent.left) {
                Node y = x.parent.parent.right;
                if (y.color == Node.RED) {
                    x.parent.color = Node.BLACK;
                    y.color = Node.BLACK;
                    x.parent.parent.color = Node.RED;
                    x = x.parent.parent;
                } else {
                    if (x == x.parent.right) {
                        x = x.parent;
                        this.leftRotate(x);
                    }
                    x.parent.color = Node.BLACK;
                    x.parent.parent.color = Node.RED;
                    this.rightRotate(x.parent.parent);
                }
            } else {
                Node y = x.parent.parent.left;
                if (y.color == Node.RED) {
                    x.parent.color = Node.BLACK;
                    y.color = Node.BLACK;
                    x.parent.parent.color = Node.RED;
                    x = x.parent.parent;
                } else {
                    if (x == x.parent.left) {
                        x = x.parent;
                        this.rightRotate(x);
                    }
                    x.parent.color = Node.BLACK;
                    x.parent.parent.color = Node.RED;
                    this.leftRotate(x.parent.parent);
                }
            }
        }
        this.root.color = Node.BLACK;
    }


    private Node root() {
        return this.root;
    }


    private Node minimum(Node node) {
        assert (node != null);
        assert (!node.isNull());
        while (!node.left.isNull()) {
            node = node.left;
        }
        return node;
    }


    private Node maximum(Node node) {
        assert (node != null);
        assert (!node.isNull());
        while (!node.right.isNull()) {
            node = node.right;
        }
        return node;
    }


    private Node successor(Node x) {
        assert (x != null);
        assert (!x.isNull());
        if (!x.right.isNull()) {
            return this.minimum(x.right);
        }
        Node y = x.parent;
        while ((!y.isNull()) && x == y.right) {
            x = y;
            y = y.parent;
        }
        return y;
    }


    private Node predecessor(Node x) {
        assert (x != null);
        assert (!x.isNull());

        if (!x.left.isNull()) {
            return this.maximum(x.left);
        }
        Node y = x.parent;
        while ((!y.isNull()) && x == y.left) {
            x = y;
            y = y.parent;
        }
        return y;
    }


    private void leftRotate(Node x) {
        Node y = x.right;
        x.right = y.left;
        if (y.left != NIL) {
            y.left.parent = x;
        }
        y.parent = x.parent;
        if (x.parent == NIL) {
            this.root = y;
        } else {
            if (x.parent.left == x) {
                x.parent.left = y;
            } else {
                x.parent.right = y;
            }
        }
        y.left = x;
        x.parent = y;

        applyUpdate(x);
        // no need to apply update on y, since it'll y is an ancestor
        // of x, and will be touched by applyUpdate().
    }


    private void rightRotate(Node x) {
        Node y = x.left;
        x.left = y.right;
        if (y.right != NIL) {
            y.right.parent = x;
        }
        y.parent = x.parent;
        if (x.parent == NIL) {
            this.root = y;
        } else {
            if (x.parent.right == x) {
                x.parent.right = y;
            } else {
                x.parent.left = y;
            }
        }
        y.right = x;
        x.parent = y;


        applyUpdate(x);
        // no need to apply update on y, since it'll y is an ancestor
        // of x, and will be touched by applyUpdate().
    }


    /**
     * Note:  Does not maintain RB constraints,  this is done post insert
     *
     * @param x
     */
    private void treeInsert(Node x) {
        Node node = this.root;
        Node y = NIL;
        while (node != NIL) {
            y = node;
            if (x.interval.getLow() <= node.interval.getLow()) {
                node = node.left;
            } else {
                node = node.right;
            }
        }
        x.parent = y;

        if (y == NIL) {
            this.root = x;
            x.left = x.right = NIL;
        } else {
            if (x.interval.getLow() <= y.interval.getLow()) {
                y.left = x;
            } else {
                y.right = x;
            }
        }

        this.applyUpdate(x);
    }


    // Applies the statistic update on the node and its ancestors.
    private void applyUpdate(Node node) {
        while (!node.isNull()) {
            this.update(node);
            node = node.parent;
        }
    }

    // Note: this method is called millions of times and is optimized for speed, or as optimized as java allows.
    private void update(Node node) {
        int nodeMax = node.left.max > node.right.max ? node.left.max : node.right.max;
        int intervalHigh = node.interval.high;
        node.max = nodeMax > intervalHigh ? nodeMax : intervalHigh;

        int nodeMin = node.left.min < node.right.min ? node.left.min : node.right.min;
        int intervalLow = node.interval.low;
        node.min = nodeMin < intervalLow ? nodeMin : intervalLow;
    }


    /**
     * Returns the number of nodes in the tree.
     */
    public int size() {
        return _size(this.root);
    }


    private int _size(Node node) {
        if (node.isNull())
            return 0;
        return 1 + _size(node.left) + _size(node.right);
    }


    private boolean allRedNodesFollowConstraints(Node node) {
        if (node.isNull())
            return true;

        if (node.color == Node.BLACK) {
            return (allRedNodesFollowConstraints(node.left) &&
                    allRedNodesFollowConstraints(node.right));
        }

        // At this point, we know we're on a RED node.
        return (node.left.color == Node.BLACK &&
                node.right.color == Node.BLACK &&
                allRedNodesFollowConstraints(node.left) &&
                allRedNodesFollowConstraints(node.right));
    }


    // Check that both ends are equally balanced in terms of black height.
    private boolean isBalancedBlackHeight(Node node) {
        if (node.isNull())
            return true;
        return (blackHeight(node.left) == blackHeight(node.right) &&
                isBalancedBlackHeight(node.left) &&
                isBalancedBlackHeight(node.right));
    }


    // The black height of a node should be left/right equal.
    private int blackHeight(Node node) {
        if (node.isNull())
            return 0;
        int leftBlackHeight = blackHeight(node.left);
        if (node.color == Node.BLACK) {
            return leftBlackHeight + 1;
        } else {
            return leftBlackHeight;
        }
    }


    /**
     * Test code: make sure that the tree has all the properties
     * defined by Red Black trees and interval trees
     * <p/>
     * o.  Root is black.
     * <p/>
     * o.  NIL is black.
     * <p/>
     * o.  Red nodes have black children.
     * <p/>
     * o.  Every path from root to leaves contains the same number of
     * black nodes.
     * <p/>
     * o.  getMax(node) is the maximum of any interval rooted at that node..
     * <p/>
     * This code is expensive, and only meant to be used for
     * assertions and testing.
     */

    public boolean isValid() {
        if (this.root.color != Node.BLACK) {
            logger.warn("root color is wrong");
            return false;
        }
        if (NIL.color != Node.BLACK) {
            logger.warn("NIL color is wrong");
            return false;
        }
        if (allRedNodesFollowConstraints(this.root) == false) {
            logger.warn("red node doesn't follow constraints");
            return false;
        }
        if (isBalancedBlackHeight(this.root) == false) {
            logger.warn("black height unbalanced");
            return false;
        }

        return hasCorrectMaxFields(this.root) &&
                hasCorrectMinFields(this.root);
    }


    private boolean hasCorrectMaxFields(Node node) {
        if (node.isNull())
            return true;
        return (getRealMax(node) == (node.max) &&
                hasCorrectMaxFields(node.left) &&
                hasCorrectMaxFields(node.right));
    }


    private boolean hasCorrectMinFields(Node node) {
        if (node.isNull())
            return true;
        return (getRealMin(node) == (node.min) &&
                hasCorrectMinFields(node.left) &&
                hasCorrectMinFields(node.right));
    }


    static class Node {

        public static boolean BLACK = false;
        public static boolean RED = true;

        Interval interval;
        int min;
        int max;
        Node left;
        Node right;

        // Color and parent are used for inserts.  If tree is immutable these are not required (no requirement
        // to store these persistently).
        boolean color;
        Node parent;


        private Node() {
            this.max = Integer.MIN_VALUE;
            this.min = Integer.MAX_VALUE;
        }

        public void store(DataOutputStream dos) throws IOException {
            dos.writeInt(interval.getLow());
            dos.writeInt(interval.getHigh());
            dos.writeInt(min);
            dos.writeInt(max);

        }

        public Node(Interval interval) {
            this();
            this.parent = NIL;
            this.left = NIL;
            this.right = NIL;
            this.interval = interval;
            this.color = RED;
        }


        static Node NIL;

        static {
            NIL = new Node();
            NIL.color = BLACK;
            NIL.parent = NIL;
            NIL.left = NIL;
            NIL.right = NIL;
        }


        public boolean isNull() {
            return this == NIL;
        }


        public String toString() {
            if (this == NIL) {
                return "nil";
            }
            /* return
                    "(" + this.interval + " " + (this.color == RED ? "RED" : "BLACK") +
                            " (" + this.left.toString() + ", " + this.right.toString() + ")";
                            */
            StringBuffer buf = new StringBuffer();
            _toString(buf);
            return buf.toString();
        }

        public void _toString(StringBuffer buf) {
            if (this == NIL) {
                buf.append("nil");
                return;
            }
            buf.append(this.interval + " -> " + this.left.interval + ", " + this.right.interval);
            buf.append("\n");
            this.left._toString(buf);
            this.right._toString(buf);
        }
    }
}

