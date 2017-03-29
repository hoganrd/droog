"""identify drug misspelling candidates using noisy channel probability model"""
#
# Copyright 2017 Robert D. Hogan
# rhogan@terminologix.com
#
# This file is part of Droog
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.

import  json
from    urllib import request
from    inspect import isfunction
from    configuration import *

# global constants
ALPHABET = 'abcdefghijklmnopqrstuvwxyz'

def generate(drug):
    """return edit distance 1 & 2 variants with noisy channel probabilities"""
    global  confusion, char, chars
    confusion = build_confusions()
    char = get_char_counts()
    chars = get_bigram_counts()
    return d1_candidates(drug), d2_candidates(drug)

def filter_variants(drug, d1, d2, limit=3, progress=False):
    """filter misspelling candidates using Google CSE"""
    candidates = d1[:int(len(d1)*D1_FILTER_PCT/100)] + \
                 d2[:int(len(d2)*D2_FILTER_PCT/100)]
    filtered = []
    filtered.append((drug, google_hits(drug)))
    for i, row in enumerate(candidates):
        filtered.append((row[0], google_hits(row[0])))
        if progress and i % 100 == 0:
            print(str(i) + '/' + str(len(candidates)))
        if limit > 0 and i > limit:     # limit=0 for all candidates
            break
    return sorted([x for x in filtered if x[1] > 0], key=lambda x: x[1], \
                  reverse=True)
        
def google_hits(word):
    """get hit count on Google for a given word (limit 100/day for free)"""
    template = 'https://www.googleapis.com/customsearch/v1?key={}&cx={}&q="{}"'
    # note: query is enclosed in quotes ("") to avoid query pre-processing
    url = template.format(CSE_KEY, CSE_CX, word)
    results = json.loads(request.urlopen(url).read().decode(), \
                         encoding='utf8')
    return int(results['queries']['request'][0]['totalResults'])

def build_confusions():
    """create confusion matrices from list of misspellings"""
    matrix = dict()
    # initialize all cells to 0
    for v_type in ['d', 'i', 's', 't']:
        matrix[v_type] = dict()
        for x in '@' + ALPHABET:
            matrix[v_type][x] = dict()
            for y in '@' + ALPHABET:
                matrix[v_type][x][y] = 0
    # process misspelling list
    with open(".\\data\\misspellings.txt", "r") as in_file:
        for line in in_file:
            misspelling, cnt, correct, v_type, x, y = line. \
                                                      replace('\n',''). \
                                                      split(',')
            cnt = int(cnt)
            matrix[v_type][x][y] = matrix[v_type][x][y] + cnt
    return matrix

def d1_candidates(word, minimum=1):
    """get candidate misspellings at edit distance = 1 with probabilities"""
    variants = insertions(word, p_insertion, minimum) + deletions(word,\
               p_deletion, minimum) + substitutions(word, p_substitution, \
               minimum) + transpositions(word, p_transposition, minimum)
    # use variants with highest probability if more than one
    v_dict = dict()
    for v in variants:
        if v[0] in v_dict:
            p = v_dict[v[0]][4]
            if p < v[4]:
                v_dict[v[0]] = v
        else:
            v_dict[v[0]] = v
    variants = v_dict.values()
    # reformat results and sort by probability (descending) and variant text (ascending)
    results = sorted([(v[0], 1, v[4], 0.0, v[1]+v[2]+v[3]) for v in variants \
                     if v[0] != word], key=lambda x: x[0])
    results = sorted([r for r in results], key=lambda x: x[2], reverse=True)
    return results

def d2_candidates(word, minimum=1):
    """get candidate misspellings at edit distance = 2 with probabilities"""
    d1s = d1_candidates(word, minimum)
    d1_dict = dict()
    for d1 in d1s:
        d1_dict[d1[0]] = True
    variants = [(d2[0], 2, d1[2], d2[2], d1[4]+'-'+d2[4]) \
                for d1 in d1s \
                for d2 in d1_candidates(d1[0], minimum) \
                if d2[0] != word and d2[0] not in d1_dict]
    # use variants with highest probability if more than one
    variants = sorted([v for v in variants], key=lambda x: (x[0], x[2]*x[3], \
                                                            x[2], x[3]))
    d2 = list()
    for i, v in enumerate(variants):
        if i+1 == len(variants) or v[0] != variants[i+1][0]:
            d2.append(v)
    # return d2 candidates sorted by joint probability (descending)
    return sorted([v for v in d2], key=lambda x: (x[2]*x[3]), reverse=True)

def insertions(word, pfunction=False, minimum=0):
    """get list of word variants with single character insertion"""
    variants = []
    for i in range(len(word) + 1):
        for j in range(len(ALPHABET)):
            y = ALPHABET[j]
            variant = word[:i] + y + word[i:]
            if i == 0:
                x = '@'
            else:
                x = word[i-1]
            if isfunction(pfunction):
                variants.append((variant, x, y, 'i', pfunction(x, y, minimum)))
            else:
                variants.append((variant, x, y, 'i'))
    return variants

def substitutions(word, pfunction=False, minimum=0):
    """get list of word variants with single character substitution"""
    variants = []
    for i in range(len(word)):
        y = word[i]
        for j in range(len(ALPHABET)):
            x = ALPHABET[j]
            variant = word[:i] + x + word[i+1:]
            if isfunction(pfunction):
                variants.append((variant, x, y, 's', pfunction(x, y, minimum)))
            else:
                variants.append((variant, x, y, 's'))
    return variants

def transpositions(word, pfunction=False, minimum=0):
    """get list of word variants with single character transposition"""
    if len(word) < 2:
        return []
    variants = []
    for i in range(2, len(word) + 1):
        x = word[i-2]
        y = word[i-1]
        variant = word[:i-2] + y + x + word[i:]
        if isfunction(pfunction):
            variants.append((variant, x, y, 't', pfunction(x, y, minimum)))
        else:
            variants.append((variant, x, y, 't'))
    return variants

def deletions(word, pfunction=False, minimum=0):
    """get list of word variants with single character deletion"""
    variants = []
    for i in range(len(word)):
        if i == 0:
            x = '@'
        else:
            x = word[i-1]
        y = word[i]
        variant = word[:i] + word[i+1:]
        if isfunction(pfunction):
            variants.append((variant, x, y, 'd', pfunction(x, y, minimum)))
        else:
            variants.append((variant, x, y, 'd'))
    return variants

def p_deletion(x, y, minimum=1, threshold=250):
    """probability of deletion of y following x"""
    global chars, confusion
    if y not in chars[x] or chars[x][y] < threshold:
        return 0.0
    else:
        return max(confusion['d'][x][y], minimum)/chars[x][y]

def p_insertion(x, y, minimum=1):
    """probability of inserting a y following x"""
    global char, confusion
    return max(confusion['i'][x][y], minimum)/char[x]

def p_substitution(x, y, minimum=1):
    """probability of substituting x for y"""
    global char, confusion
    return max(confusion['s'][x][y], minimum)/char[y]

def p_transposition(x, y, minimum=1, threshold=250):
    """probability of replacing xy with yx"""
    global chars, confusion
    if y not in chars[x] or chars[x][y] < threshold:
        return 0.0
    else:
        return max(confusion['t'][x][y], minimum)/chars[x][y]

def get_char_counts():
    """return dictionary of characters counts"""
    counts = read_delimited_list('.\\data', 'char_counts.txt')
    result = dict()
    for element in counts:
        result[element[0]] = int(element[1])
    return result

def get_bigram_counts():
    """return dictionary of bigram counts"""
    counts = read_delimited_list('.\\data', 'bigram_counts.txt')
    result = dict()
    for element in counts:
        bigram = element[0]
        count = element[1]
        l1 = bigram[0]
        l2 = bigram[1]
        if l1 not in result:
            result[l1] = dict()
        result[l1][l2] = int(count)
    return result

def read_delimited_list(path, file_name):
    """read comma separated list of character or bigram counts"""
    result = []
    with open(path + '\\' + file_name, 'r') as file:
        for line in file:
            result.append(line.strip().split(','))
    return result
