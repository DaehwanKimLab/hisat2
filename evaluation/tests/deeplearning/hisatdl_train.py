#!/usr/bin/env python
#
# Copyright 2017, Daehwan Kim <infphilo@gmail.com>
#
# This file is part of HISAT 2.
#
# HISAT 2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT 2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.
#


import os, sys, subprocess, re
import tensorflow as tf
import numpy as np
from argparse import ArgumentParser, FileType


"""
# The following example taken from https://www.tensorflow.org/get_started/get_started
"""
def tf_example():
    # Model parameters
    W = tf.Variable([.3], dtype=tf.float32)
    b = tf.Variable([-.3], dtype=tf.float32)
    # Model input and output
    x = tf.placeholder(tf.float32)
    linear_model = W*x + b
    y = tf.placeholder(tf.float32)

    # loss
    loss = tf.reduce_sum(tf.square(linear_model - y)) # sum of the squares
    # optimizer
    optimizer = tf.train.GradientDescentOptimizer(0.01)
    train = optimizer.minimize(loss)

    # training data
    x_train = [1, 2, 3, 4]
    y_train = [0, -1, -2, -3]
    # training loop
    init = tf.global_variables_initializer()
    sess = tf.Session()
    sess.run(init) # reset values to wrong
    for i in range(1000):
      sess.run(train, {x: x_train, y: y_train})

    # evaluate training accuracy
    curr_W, curr_b, curr_loss = sess.run([W, b, loss], {x: x_train, y: y_train})
    print("W: %s b: %s loss: %s"%(curr_W, curr_b, curr_loss))


"""
"""
def read_seqs(fname):
    feature_len = None
    features = []
    for line in open(fname):
        line = line.strip()
        feature = []
        for num in line.split(' '):
            if num == '0':
                feature.append(0.0)
            else:
                feature.append(1.0)
        if feature_len != None:
            assert feature_len == len(feature)
        else:
            feature_len = len(feature)
        features.append(feature)
    return features


"""
"""
def read_labels(fname):
    label_len = None
    labels = []
    for line in open(fname):
        line = line.strip()
        label = []
        for num in line.split(' '):
            label.append(float(num))
        if label_len != None:
            assert label_len == len(label)
        else:
            label_len = len(label)
        labels.append(label)
    return labels

def weight_variable(shape):
  initial = tf.truncated_normal(shape, stddev=0.1)
  return tf.Variable(initial)

def bias_variable(shape):
  initial = tf.constant(0.1, shape=shape)
  return tf.Variable(initial)

def conv1d(x, W):
  return tf.nn.conv1d(x, W, stride=1, padding='SAME')

def max_pool_2(x):
  return tf.nn.max_pool(x, ksize=[1, 2, 1],
                        strides=[1, 2, 1], padding='SAME')


"""
"""
def train_sequences_CNN(base_fname,
                        locus_list,
                        verbose):
    # filenames = ["label.train", "seq.train"]
    # dataset = tf.contrib.data.TextLineDataset(filenames)
    # print dataset

    batch_size = 1000
    features_array, labels_array = read_seqs("seq.train"), read_labels("label.train")
    feature_len, label_len = len(features_array[0]), len(labels_array[0])
    assert len(features_array) > 0 and len(features_array) == len(labels_array)

    x = tf.placeholder(tf.float32, [None, feature_len])
    x_reshape = tf.reshape(x, [-1, feature_len, 1])
    y = tf.placeholder(tf.float32, [None, label_len])

    W_conv1 = weight_variable([12, 1, 8])
    b_conv1 = bias_variable([8])
    h_conv1 = tf.nn.relu(conv1d(x_reshape, W_conv1) + b_conv1)
    # h_pool1 = max_pool_2(h_conv1)

    W_conv2 = weight_variable([12, 8, 16])
    b_conv2 = bias_variable([16])
    h_conv2 = tf.nn.relu(conv1d(h_conv1, W_conv2) + b_conv2)

    W_conv3 = weight_variable([12, 16, 32])
    b_conv3 = bias_variable([32])
    h_conv3 = tf.nn.relu(conv1d(h_conv2, W_conv3) + b_conv3)

    W_conv4 = weight_variable([12, 32, 64])
    b_conv4 = bias_variable([64])
    h_conv4 = tf.nn.relu(conv1d(h_conv3, W_conv4) + b_conv4)

    W_conv5 = weight_variable([12, 64, 64])
    b_conv5 = bias_variable([64])
    h_conv5 = tf.nn.relu(conv1d(h_conv4, W_conv5) + b_conv5)

    W_fc1 = weight_variable([40 * 64, 64])
    b_fc1 = bias_variable([64])

    h_conv2_flat = tf.reshape(h_conv5, [-1, 40 * 64])
    h_fc1 = tf.nn.relu(tf.matmul(h_conv2_flat, W_fc1) + b_fc1)

    W_fc2 = weight_variable([64, 64])
    b_fc2 = bias_variable([64])
    h_fc2 = tf.nn.relu(tf.matmul(h_fc1, W_fc2) + b_fc2)

    keep_prob = tf.placeholder(tf.float32)
    keep_prob = 1.0
    h_fc1_drop = tf.nn.dropout(h_fc1, keep_prob)

    W_fc3 = weight_variable([64, 10])
    b_fc3 = bias_variable([10])

    y_ = tf.matmul(h_fc2, W_fc3) + b_fc3
    ysoftmax = tf.nn.softmax(y_)
    cross_entropy = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(labels=y, logits=y_))
    optimizer = tf.train.AdamOptimizer(1e-4)
    train = optimizer.minimize(cross_entropy)

    sess = tf.Session()
    init = tf.global_variables_initializer()
    sess.run(init)
    for i in range(100000):
        for b in range(len(features_array) / batch_size + 1):
            features, labels = features_array[b*batch_size:(b+1)*batch_size], labels_array[b*batch_size:(b+1)*batch_size]
            if len(features) <= 0:
                break
            sess.run(train, {x: features, y: labels})
        if i % 500 == 0:
            print "DK:", i
            # print "\t", sess.run([cross_entropy], {x: features_array, y: labels_array})

            # evaluate training accuracy
            # correct_prediction = tf.equal(tf.argmax(y,1), tf.argmax(y_,1))
            # accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
            # print "\t", sess.run(accuracy, {x: features_array, y: labels_array})
    
            ucorrect_count, utotal_count = 0, 0
            mcorrect_count, mtotal_count = 0, 0
            for b in range(len(features_array) / batch_size + 1):
                features, labels = features_array[b*batch_size:(b+1)*batch_size], labels_array[b*batch_size:(b+1)*batch_size]
                if len(features) <= 0:
                    break
                predicted_labels = sess.run(ysoftmax, {x: features, y: labels})
                for l in range(len(labels)):
                    label, predicted_label = labels[l], predicted_labels[l]
                    if max(label) == 1.0:
                        label_max = np.argmax(label)
                        predicted_label_max = np.argmax(predicted_label)
                        if label_max == predicted_label_max:
                            ucorrect_count += 1
                        utotal_count += 1
                    else:
                        same = True
                        for i_ in range(len(label)):
                            if label[i_] > 0.0:
                                # if abs(label[i_] - predicted_label[i_]) / label[i_] > 0.7:
                                if predicted_label[i_] < 0.02:
                                    same = False
                        if same:
                            mcorrect_count += 1
                        mtotal_count += 1

            print >> sys.stderr, "\tunique: %d / %d (%.2f), multi: %d / %d (%.2f)" % \
                (ucorrect_count, utotal_count, float(ucorrect_count) / utotal_count * 100, mcorrect_count, mtotal_count, float(mcorrect_count) / mtotal_count * 100)


    
"""
"""
def train_sequences_W(base_fname,
                      locus_list,
                      verbose):
    # filenames = ["label.train", "seq.train"]
    # dataset = tf.data.TextLineDataset(filenames)
    # print dataset

    features, labels = read_seqs("seq.train"), read_labels("label.train")
    feature_len, label_len = len(features[0]), len(labels[0])
    assert len(features) > 0 and len(features) == len(labels)

    x = tf.placeholder(tf.float32, [None, feature_len])
    y = tf.placeholder(tf.float32, [None, label_len])
    W = tf.Variable(tf.zeros([feature_len, label_len]))
    b = tf.Variable(tf.zeros([label_len]))
    y_ = tf.nn.softmax(tf.matmul(x,W) + b)

    cross_entropy = tf.reduce_mean(-tf.reduce_sum(y * tf.log(y_), reduction_indices=[1]))
    optimizer = tf.train.GradientDescentOptimizer(0.01)
    train = optimizer.minimize(cross_entropy)
    
    sess = tf.Session()
    init = tf.global_variables_initializer()
    sess.run(init)
    for i in range(100):
        sess.run(train, {x: features, y: labels})
        if i % 10 == 0:
            print "DK:", i
            print "\t", sess.run([cross_entropy], {x: features, y: labels})

    # evaluate training accuracy
    curr_W, curr_b, curr_loss = sess.run([W, b, cross_entropy], {x: features, y: labels})
    print("W: %s b: %s loss: %s"%(curr_W, curr_b, curr_loss))

    correct_prediction = tf.equal(tf.argmax(y,1), tf.argmax(y_,1))
    accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
    print sess.run(accuracy, {x: features, y: labels})

    
        
"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description="Generate train sequences")
    parser.add_argument("-b", "--base",
                        dest="base_fname",
                        type=str,
                        default="hla",
                        help="base filename for backbone sequence, variants, and linking info (Default: hla)")
    parser.add_argument("--locus-list",
                        dest="locus_list",
                        type=str,
                        default="",
                        help="A comma-separated list of gene names (default: empty, all genes)")
    parser.add_argument("-v", "--verbose",
                        dest="verbose",
                        action="store_true",
                        help="also print some statistics to stderr")

    args = parser.parse_args()
    if args.locus_list == "":
        locus_list = []
    else:
        locus_list = args.locus_list.split(',')
             
    if args.base_fname.find('/') != -1:
        elems = args.base_fname.split('/')
        base_fname = elems[-1]
        base_dname = '/'.join(elems[:-1])
    else:
        base_fname = args.base_fname
        base_dname = ""

    train_sequences_CNN(base_fname,
                        locus_list,
                        args.verbose)

