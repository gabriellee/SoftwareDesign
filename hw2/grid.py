# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 18:21:43 2014

@author: gabrielle
"""

def grid(x,y):
    i=1
    j=1
    while i<=y:
        print ('+' + '-'*4)*x + '+'
        
        print ("|" +' '*4)*x + '|'
        print ("|" +' '*4)*x + '|'
        print ("|" +' '*4)*x + '|'
        print ("|" +' '*4)*x + '|'
        i=i+1
    print ('+' + '-'*4)*x + '+'        
        
            
