#!/usr/bin/env python

import sys, os, subprocess, re 
import inspect, random 
import math 
import glob
import urllib2
import string

if __name__ == '__main__':

    result_list = {}
    alleles = []
    calls = []
    
    for line in open('ABO_syn_results_paired.txt', 'r'):
        line = line.replace('***', '').strip()
        if 'Test' in line:
            if alleles and calls:
                for allele in alleles:
                    for call in calls:
                        if call not in result_list[allele][2]:
                            result_list[allele][2].update({ call : 1 })
                        else:
                            result_list[allele][2][call] += 1
            alleles = []
            calls = []

        if line.startswith('ABO'):
            allele = line.split()[0]
            if allele not in result_list:
                result_list.update({ allele : [1,0,{}] })
            else:
                result_list[allele][0] += 1
            alleles.append(allele)

        if 'abundance' in line:
            if line.startswith('1 ') or line.startswith('2 '):
                allele = line.split()[2]
                if allele in alleles:
                    result_list[allele][1] += 1
                    alleles.remove(allele)
                else:
                    calls.append(allele)

    
    for allele, counts in result_list.items():
        print '%s\tTotal: %d\tCorrect: %d\tPercent: %d\tCalled as: %s' % (allele, counts[0], counts[1], (counts[1]*100)/counts[0], counts[2])         
