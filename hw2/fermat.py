# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 19:30:16 2014

@author: gabrielle
"""
print "Please enter the value for a"
a = int(raw_input())
print "Please enter the value for b"
b = int(raw_input())
print "Please enter the value for c"
c = int(raw_input())
print "Please enter the value for n"
n = int(raw_input())

def check_fermat(a,b,c,n):
    if n>2 and a**n + b**n == c**n:
        print "Holy smokes, Fermat was wrong!"
    else:
        print "No, that doesn't work"
        
        
check_fermat(a,b,c,n)