#!/usr/bin/env python

import sys
import subprocess
import os
import re
import pickle
import datetime
import time
import signal

executable = '../cave'

vars = 'abcdefghijklmnopqrstuvwxyz'

TIMEOUT = 500           #secs to wait for result
TIME_INTERVAL = 0.1     #secs between polls
SAVE_FREQ = TIMEOUT    #save every this-many seconds (we can't save a partial result though)
                       # so SAVE_FREQ = TIMEOUT will save after each test

def timeout_command(command, timeout):
    """call shell-command and either return its output or kill it
    if it doesn't normally exit within timeout seconds and return None"""
    start = datetime.datetime.now()
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    t = 0.0
    while (process.poll() is None):
        t = t + 0.1
        time.sleep(0.1)
        now = datetime.datetime.now()
        if ((now - start).seconds > timeout):
            os.kill(process.pid, signal.SIGKILL)
            os.waitpid(-1, os.WNOHANG)
            return None
            
    return (t, process.stdout.read())

def indegree_at_least_formula(k):
    formula = "'A" + vars[0] + ' '
    #quantifiers
    for i in range(k):
        formula += 'E'
        formula += vars[i+1]
        formula += ' '
    formula += '('
    for i in range(k):
        formula += '('
        formula += vars[i+1]
        formula += '->'
        formula += vars[0]
        formula += ')&'
    for i in range(k):
        for j in range(i):
            if (i != j):
                formula += '~('
                formula += vars[i+1]
                formula += '=='
                formula += vars[j+1]
                formula += ')&'
    formula = formula[:-1] + ')' + "'"
    return formula;


def k_cycle_formula(k):
    formula = "'"
    #quantifiers
    for i in range(k):
        formula += 'E'
        formula += vars[i]
        formula += ' '
    formula += '('
    for i in range(k):
        formula += '('
        formula += vars[i]
        formula += '->'
        formula += vars[(i+1)%k]
        formula += ')&'
    for i in range(k):
        if (i > 0) and (k % i) == 0:
            formula += '~('
            formula += vars[0]
            formula += '=='
            formula += vars[i]
            formula += ')&'
    formula = formula[:-1] + ')' + "'"
    return formula;

def check_formula(formula, eca, omega):
    try:
        args = [executable, '-e', '%d' % eca, '-f ', formula ]
        if not omega:
            args.append(' -Z')
        #p = Popen(' '.join(args), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output = timeout_command(' '.join(args,), TIMEOUT) 
        if(output is None):
            return ('TIMEOUT@%d'% TIMEOUT, None) 
        sts = (output[1].find('true') != -1)
        if output[1].find('Error parsing formula') != -1:
            print "Bad things happened", -retcode
            return None
        else:
            return (output[0], sts)
    except OSError, e:
        print "Execution failed:", e
        return None

def load(savefile):
    fileobj = open(savefile, 'r')
    ret = {}
    for line in fileobj:
        parts = line.split(' ')
        if(len(parts) > 3):
            key = (int(parts[0]), int(parts[1]))
            val = (parts[2], parts[3])
            ret[key] = val
    fileobj.close()       
    print 'Loaded:'
    print ret
    return ret
    
def saved(d, savefile, timestamp):
    if(time.time() - timestamp > SAVE_FREQ):
        fileobj = open(savefile, 'w')
        for (k, eca) in d:
            fileobj.write('%d %d ' % (k, eca))
            value = d[(k, eca)]
            if type(value[0]).__name__=='float':
                fileobj.write( '%3.1f %s\n' % value )
            else:
                fileobj.write( '%s %s\n' % value )            
        
        fileobj.close()
        return time.time()
    else:
        return timestamp
        
def main():
    MAX_K = 6
    OMEGA = True
    
    savefile = './indegree'
    
    K = map(lambda x : x + 2, range(MAX_K-1))
    
    for arg in sys.argv:
        if arg.find('zeta') != -1:
            OMEGA = False   

    k_eca_pairs = [(k,eca) for k in K for eca in range(256)]
    
    d = {}
    if(os.path.isfile(savefile)):
        d = load(savefile)
        
    timestamp = time.time()
    for (k, eca) in k_eca_pairs:
        # We need to check this pair if: 
        # 1) It has not yet been tested
        check = (k, eca) not in d
        # or 2) it has been tested and timed out with a timeout that was
        # less than our current value.
        if(not check):
            if (d[(k,eca)][1] == 'None\n'):
                old_timeout = int(d[(k,eca)][0].split('@')[1])
                if(old_timeout < TIMEOUT):
                    check = True
        # check the pair if one of the above conditions were met                    
        if (check): 
            formula = indegree_at_least_formula(k)
            print 'Checking %d %d' % (k, eca)
            result = check_formula(formula, eca, OMEGA)
            if(result != None):
                d[(k,eca)] = result
                print result
                timestamp = saved(d, savefile, timestamp)
                
main()
