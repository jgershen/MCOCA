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

TIMEOUT = 120           #secs to wait for result
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

def check_k_cycle(k, eca, omega):
    try:
        args = [executable, '-e', '%d' % eca, '-f ', k_cycle_formula(k) ]
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
        
def output(data):
    print data
    outf = open('./output', 'a')
    outf.write(data)
    outf.flush()
    outf.close()
    
def main():
    MAX_K = 16
    OMEGA = True
    
    savefile = './cycles'
    
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
        if (k, eca) not in d or d[(k,eca)][1] is None: #recheck if a previous timeout
            result = check_k_cycle(k, eca, OMEGA)
            if(result != None):
                d[(k,eca)] = result
                if type(result[0]).__name__=='float':
                    output( '%d %d %3.1f %s\n' % (k, eca, result[0], result[1]) )
                else:
                    output( '%d %d %s %s\n' % (k, eca, result[0], result[1]) )
                timestamp = saved(d, savefile, timestamp)
    
    print d
    
    '''
    vals = [[0 for x in range(256)] for y in range(MAX_K-1)]

    for k in range(MAX_K):
        if(k >= 1):
            for eca in range(256):
                if(check_k_cycle(k, eca, OMEGA)):
                    vals[k-1][eca] = 1
                    
    print vals
    '''
main()
