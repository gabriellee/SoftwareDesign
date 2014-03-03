# -*- coding: utf-8 -*-
"""
Created on Sat Mar  1 17:53:00 2014

@author: gabrielle
"""
def get_lingual_distance(book):
    #book is a translation dictionary
    #word_pair_list contains a list of lists for each word defined in the dictionary
    word_pair_list = parse_word_pairs(book)
    distances = list()
    num_subst = 0
    num_ins = 0
    num_del = 0
    for i in range(len(word_pair_list)):
        for j in range(len(word_pair_list[i])):
            output = levenshtein_fun(word_pair_list[i][0],word_pair_list[i][j])
            distances.append(output[0])
            num_subst += output[1]
            num_ins += output[2]
            num_del += output[3]
    avg_dist = 
    return 
            
def parse_word_pairs():
def levenshtein_fun(s1, s2):
    num_subst = 0
    num_ins = 0
    num_del = 0
    #s2 must be the shorter string
    if len(s1) < len(s2):
        temp = s1;
        s1 = s2;
        s2 = temp;
    if len(s2) == 0:
        return len(s1)
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
            deletions = current_row[j] +1
            #add one to the edit distance between the prior c2 and the prior c1
            substitutions = prior_row[j] + 1
            current_row.append(min(insertions, deletions, substitutions))
            if min(insertions, deletions, substitutions) == substitutions:
                num_subst+=1
            elif min(insertions, deletions, substitutions) == insertions:
                num_ins+=1
            else:
                num_del+=1
        prior_row = current_row
    return prior_row, num_subst, num_ins, num_del
            
        