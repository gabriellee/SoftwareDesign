# -*- coding: utf-8 -*-
"""
Created on Sat Mar  1 17:53:00 2014

@author: gabrielle
"""
import re
import pickle
g = open('German.txt','r')
german = g.read()
g.close()
f = open('LengthbySize.pickle','w')
def parse_word_pairs(book):
    word_list = re.split('\r\n\r\n',book)
    word_list.remove('\xef\xbb\xbf')
    pairs = []
    for i in range(len(word_list)):
        s = re.split('  +',word_list[i])
        pairs.append((s[0].split(),s[1].split()))
    return pairs

def levenshtein_fun(s1, s2):
    num_subst = 0
    num_ins = 0
    num_del = 0
    #s2 must be the shorter string
    if len(s1) < len(s2):
        return levenshtein_fun(s2, s1)
    if len(s2) == 0:
        return len(s1), 0, len(s1), 0
    #the inital row is a comparison between the first character of s1 and each character of s2
    prior_row = range(len(s2)+1)
    #loop through each character in s1(the shorter string) and compare it to consecutive characters in s2
    #construct the current row by adding the min edit-distance for each character of s2
    #c1 is the character in string 1
    #c2 is the character in string 2
    for i,c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            #the value of a c1 to c2 comparison is dependent upon the comparison of the prior c1 to the same c2 (edit distance is cumulative)
            #insertions, deletions, and substitutions all refer to edit-distances between s1[:i] and s2[:j]
            #transforming a string by insertion has an edit distance of 1
            insertions = prior_row[j + 1] + 1
            #add 1 to the edit distance between the prior c2 and the same c1
            deletions = current_row[j] + 1
            #add one to the edit distance between the prior c2 and the prior c1
            substitutions = prior_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        prior_row = current_row
        weighted_dist = float(prior_row[-1])/float(len(s1))
    return weighted_dist

def get_lingual_distance(book):
    #book is a translation dictionary
    #word_pair_list contains a list of lists for each word defined in the dictionary
    #the output is a long string that contains the levenshtein distance, two spaces, and then the length of the string
    #This is then put in matlab and plotted
    word_pair_list = parse_word_pairs(book)
    distances = list()
    num_subst = 0
    num_ins = 0
    num_del = 0
    out = []
    for i in range(len(word_pair_list)):
        dist =  []
        leg = []
        german = word_pair_list[i][0]
        english = word_pair_list[i][1]
        for j in range(len(german)):
            for k in range(len(english)):
                dist.append(levenshtein_fun(german[j],english[k]))
                leg.append(max([len(german[j]) , len(english[k])]))
        out.append([str(min(dist)),str(max(leg))])
        outstring = listtostr(out)
        pickle.dump("apples",f)
        f.close()

def listtostr(s):
    string = ""
    for i in range(len(s)):
        for j in range(len(s[i])):
            string = string + s[i][j] + '  '
        string = string + '\n'
    return string 





            