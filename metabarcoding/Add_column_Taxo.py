#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
#author          :Annie Lebreton
#date            :01/03/2020
#version         :1
#python_version  :3.6
"""

import argparse

def parseoptions( ):
    """ Docstring
    .... """

    parser = argparse.ArgumentParser( description="aim to split a taxonomy column with the format k__XXX;p__XXX... into multiples columns" )
    parser.add_argument( '-i',  '--infile', default="ITS_Ailaoshan_abundance.tsv",  help="taxonomy in second column (1rst in python)")

    global ARGS         # Update the global ARGS variable
    ARGS = parser.parse_args()

##----- INPUT ---------

##----- OUTPUT --------


def main():
        dico={}
        lst=["Kingdom","Phylum", "Class","Order","Family","Genus","Species"]
        lst2=["k__","p__","c__","o__","f__","g__","s__"]

        for i in lst:
            dico[i] = "Unknown"

        c = -1
        parseoptions()
        with open(ARGS.infile, 'r') as ef:
            lines = ef.readlines()
            for line in lines:
                line = line.strip()
                lS = line.split('\t')
                if c == -1:
                    while c < len(lst) - 1:
                        c = c + 1
                        dico[lst2[c]] = lst[c]

                    l = lS + [dico["k__"], dico["p__"], dico["c__"], dico["o__"], dico["f__"], dico["g__"], dico["s__"]]
                    print("\t".join(l))
                    continue

                taxo = lS[1]
                lS2 = taxo.split(";")
                c2 = 0
                l = lS
                while c2 < len(lS2):
                    tmp = lS2[c2].split("__")
                    val = tmp[-1]
                    l = l + [val]
                    c2 = c2 + 1

                print("\t".join(l))


main()
