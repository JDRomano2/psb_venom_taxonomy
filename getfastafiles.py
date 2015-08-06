#!/usr/bin/env python

import urllib2
import pdb
import re
import subprocess
import shlex # Format command line args
from tempfile import NamedTemporaryFile
import time
import json
import csv
from itertools import product

pfams = [
    "PF00008",
#    "PF00033",
#    "PF00042",
##    "PF00067",
    "PF00069",
    "PF00119",
#    "PF00207",
    "PF00243",
##    "PF02931",
    "PF06444",
    "PF00115"
]


def fasta_reader(fasta_text):
    
    fasta_dict = {}
    split_text = fasta_text.split('>')
    split_text.pop(0) # First element is empty
    for f in split_text:
        curdict = {}
        ls = f.splitlines()
        header = ls[0]
        seq = "".join(ls[1:])
        try:
            #FOR SPECIES:
            #species = re.search("OS=('?[A-z\-]+'?\s\(?'?[A-z]+)", header).group(1)
            #FOR GENUS:
            species = re.search("OS=('?[A-z\-]+)'?\s\(?'?[A-z]+", header).group(1)
     
        except Exception:
            print header
            continue # CONTINUE IS IMPORTANT!! `pass` is basically a useless function...
        try:
            curdict['accession'] = re.match("^[a-z]+\|([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})\|", header).group(1)
        except Exception:
            #pdb.set_trace()
            print header
            continue
        curdict['nameline'] = header
        curdict['seq'] = seq
        
        if species in fasta_dict:
            fasta_dict[species].append(curdict)
        else:
            fasta_dict[species] = [curdict]
        # clean up
        species = None
        curdict = None
    return fasta_dict

def get_best_seqs(all_fasta_dicts, intersection):
    
    # Find top scoring "Homo sapiens" sequence
    
    hsap = [ {k: v['Xenopeltis'][0]} for k,v in all_fasta_dicts.iteritems() ] # <== NOTE!!! This 'roots' the analysis. May want to change at some point to a more optimal taxon
    # hsap = []
    # for x in all_fasta_dicts:
    #     pdb.set_trace()
    #     hsap.append(x['Homo sapiens'][0])
    hsap = dict(zip(all_fasta_dicts.keys(), hsap))

    # For each "best" human sequence, find best alignment for each species in `intersection`
    pfam_data = {}
    for k,v in hsap.iteritems():
        best_seqs = {'hsap': v[k]}
        for i in intersection:
            if len(all_fasta_dicts[k][i]) == 1:
                best_seq = all_fasta_dicts[k][i][0]['accession']
            else:
                best_seq = doblastp(v[k], all_fasta_dicts[k][i])
            best_seq_dict = [ x for x in all_fasta_dicts[k][i] if x['accession'] == best_seq ][0]
            best_seqs[i] = best_seq_dict
        pfam_data[k] = best_seqs
    #pdb.set_trace()
    #print ""
    return pfam_data
            

def doblastp(query, db_seqs):
    # Build blast database
    
    blastdbstr = ""
    for d in db_seqs:
        blastdbstr += ">{0}\n".format(d['nameline'])
        blastdbstr += "{0}\n".format(d['seq'])
    subject = NamedTemporaryFile()
    nm = subject.name
    subject.write(blastdbstr)
    subject.seek(0)
    # command_line = "makeblastdb -in {0} -out {1}.fa -parse_seqids -dbtype prot".format(subject.name, subject.name)
    # args = shlex.split(command_line)
    args = ["makeblastdb", "-in", subject.name, "-out", "{0}.fa".format(subject.name), "-parse_seqids", "-dbtype", "prot", "-logfile", "/dev/null"]
    
    subprocess.Popen(args)
    time.sleep(0.1) #apparently I need to wait for the makeblastdb process to complete!
        
    querystr = ">{0}\n{1}\n".format(query['nameline'], query['seq'])
    query = NamedTemporaryFile()
    query.write(querystr)
    query.seek(0)

    args = None
    while True:
        try:
            command_line = "blastp -query {0} -db {1}.fa -outfmt 6".format(query.name, subject.name)
            args = shlex.split(command_line)
            result = subprocess.check_output(args)
            if result == "":
                # It's possible that no hits are found, for whatever reason. If that's the case, just return the first sequence as the best
                acc = db_seqs[0]['accession']
                return acc
            break
        except Exception:
            pdb.set_trace()
            print ""

        
    # Parse result
    # It seems the first line is highest scoring (and best)
    try:
        top_res = result.split("\n")[0]
        acc = re.match("^\w+\|\w+\|\w+\\t\w+\|(\w+)\|.+$", top_res).group(1)
    except Exception:
        print "Error parsing BLASTp results"
        pdb.set_trace()
        print ""


    subject.close()
    query.close()
    # clean up tmp directory
    command_line = "rm {0}*".format(nm)
    args = None
    args = shlex.split(command_line)
    subprocess.Popen(command_line, shell=True) # No wildcard expansion unless through shell (or see http://stackoverflow.com/questions/7156892/wildcard-not-working-in-subprocess-call-using-shlex)
    
    return acc

# def filter(l):
#     # load CSV file:
#     with open("reptiles.csv", 'rb') as reptcsv:
#         reptilereader = csv.reader(reptcsv, delimiter=',')
#         reptiles = []
#         for row in reptilereader:
#             reptiles.append(row)
#         headers = reptiles[0]
#         del(reptiles[0])
#         allrepts_list = [ {x[0]: x[2]} for x in reptiles ]
#         allrepts = {}
#         for arl in allrepts_list:
#             #pdb.set_trace()
#             allrepts[dict.keys(arl)[0]] = dict.values(arl)[0]


#     # first see how many of our species are reptiles
#     numrepts = 0
#     for r in l:
#         if r in dict.keys(allrepts):
#             print r
            
#     # scan through l and see if match in allrepts has value containing 'squamata'
#     #indices = []
#     #for i in l:
#     #    if "squamata" in allrepts[i]
        
        
#     pdb.set_trace()
#     print ""

def main():
    all_fasta_dicts = {}
    for pfam in pfams:
        # Load FASTA file
        
        with open ("raw_fasta/reptiles/{0}.fasta".format(pfam.lower()), "r") as fastafile:
            raw_fasta = fastafile.read()

        #pdb.set_trace()
        # Parse FASTA file (store FASTA sequences by species)
        fasta_dict = fasta_reader(raw_fasta)

        

        all_fasta_dicts[pfam] = fasta_dict

    # Analysis of species spanning all:
    intersection = []

    #pdb.set_trace()
    print ""
    for k in all_fasta_dicts:
        
        if len(intersection) == 0:
            intersection = set(all_fasta_dicts[k].keys())
            print "THIS SHOULD ONLY PRINT ONCE!!!"
        else:
            intersection = intersection & set(all_fasta_dicts[k].keys())
            if (len(intersection) == 0):
                raise Exception

                    
    
    intersection = list(intersection)

    ## Constrain list to only squamates
    #intersection = filter(intersection)

    #pdb.set_trace()
    #print ""

    # Construct dict of sequences to hold for analysis
    best_seqs = get_best_seqs(all_fasta_dicts, intersection)

    # Align sequences
    #aligned_seqs = align_sequences(best_seqs)

    print "Got to the end!!!!"
    pdb.set_trace()
    print ""
    
    with open('protdata.json', 'w') as fp:
        json.dump(best_seqs, fp)
            

if __name__ == "__main__":
    main()
