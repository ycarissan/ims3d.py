#!/usr/bin/python3
# encoding utf-8

import os
import argparse
import subprocess
import logging
import random
import rdkit.Chem

class MolecularGraph():

    def __init__(self,geomfile):
        mol_filename = generate_mol(geomfile, None)
        self.molecule = rdkit.Chem.MolFromMolFile(mol_filename)

    def getCycles(self):
        return self.molecule.GetRingInfo().AtomRings()

    def getEdges(self):
        return self.molecule.GetBonds()

def generate_mol(geomfile, logger):
    """ Genere le fichier de descritpion de la geometrie en mol"""
    mol_filename = "tmpfile_{:05d}.mol".format(int(random.uniform(0, 99999)))
    if (logger):
        logger.info("Temporary mol file : {}".format(mol_filename))
    proc = subprocess.Popen(["obabel",
                             "-ixyz",
                             geomfile,
                             "-omol",
                             "-O",
                             mol_filename],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,
                            universal_newlines=True)
    stdout, stderr = proc.communicate()
    if (logger):
        # stderr is None because of stderr=subprocess.STDOUT in the Popen call
        logger.info(stdout)
    return mol_filename

def main():
    print("Cycle detection library")


if __name__ == "__main__":
    main()
