#!/usr/bin/env python

# Author: Carly Hamilton carlyham@uoregon.edu

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the UO Bioinformatics and Genomics Program, 2024 coursework.
It includes functions to evaluate phred scores, calculate the median of a sorted list,
create a one line fatsa file, and more. Last updated July 19th, 2024'''

__version__ = "0.5"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = set('ATGCNatcgn')
RNA_bases = set('AUGCNaucgn')

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return(ord(letter)-33)

def qual_score(phred_score: str) -> float:
    '''Calculates an average quality score from a phred string'''
    sum_scores = 0
    for letter in phred_score:
        score = convert_phred(letter) 
        sum_scores += score
    return (sum_scores/len(phred_score))  

def validate_base_seq(seq: str, RNAflag: bool=False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNA_bases if RNAflag else DNA_bases)

def gc_content(seq: str):
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1. Case insensitive'''
    assert validate_base_seq(seq) == True, "String contains invalid characters - are you sure you used a DNA sequence?"
    seq = seq.upper()
    Gs = seq.count("G")
    Cs = seq.count("C")
    return (Gs+Cs)/len(seq)

def calc_median(lst: list) -> float:
    '''Given a sorted list, returns the median value of the list'''
    n = len(lst)
    if len(lst)%2==0:
        mid1 = int((n-1)/2) 
        mid2 = int(((n-1)/2)+1)
        med = (lst[mid1] + lst[mid2])/2
    else:
        index = int(n/2)
        med = lst[index]
    return med

def oneline_fasta(filein: str, fileout:str):
    '''Turns multiple sequence lines into a single line 
    per header in a fasta file'''
    first_line = True
    with open(filein, "r") as fin:
        with open(fileout, "w") as fout:
            for line in fin:
                line = line.strip("\n")
                if line[0] == ">":
                    if first_line == True:
                        fout.write(f"{line}\n")
                        first_line = False
                    else:
                        fout.write(f"\n{line}\n")
                else:
                    fout.write(line)

if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")

    #qual score assertions
    assert qual_score("A") == 32.0, "wrong average phred score for 'A'"
    assert qual_score("AC") == 33.0, "wrong average phred score for 'AC'"
    assert qual_score("@@##") == 16.5, "wrong average phred score for '@@##'"
    assert qual_score("EEEEAAA!") == 30.0, "wrong average phred score for 'EEEEAAA!'"
    assert qual_score("$") == 3.0, "wrong average phred score for '$'"
    print("Quality score function is working!")

    #gc content assertions
    assert gc_content("GGCGTAAT") == 0.5, "wrong GC content calculated for DNA"
    print("gc_content function is working!")

    #Assertion tests for validate base seq
    assert validate_base_seq("CTAAGGACT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("UAGCAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonRNA"
    print("Passed DNA and RNA tests")

    #Assertion tests for calc_median
    assert calc_median([3,5,10]) == 5, "calc_median function does not work for odd length list"
    assert calc_median([4,5]) == 4.5, "calc_median function does not work for even length list"
    print("Median successfully calculated")
    



    
