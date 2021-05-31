'''
@author: Oliver Unke
'''
import copy
from openbabel import openbabel as ob

obConversion = ob.OBConversion()

def traverse(idx,graph,max_size,visited=None):
    if max_size < 1:
        return
    if visited is None:
        visited = set()
    else:
        visited = copy.deepcopy(visited)
    visited.add(idx)
    yield visited
    if max_size > 1:
        for i in graph[idx]:
            if i not in visited:
                for subgraph in traverse(i,graph,max_size-1,visited):
                    yield subgraph

def bottom_up_amons_of(mol, max_size, already_generated=set()):
    """
    Generates all unique amons of mol by growing graphs up to a max_size around each atom
    already_generated: A set of smiles strings of amons that have already been generated
    """
    if mol.NumHvyAtoms() < 1:
        return

    obConversion.SetInAndOutFormats("smi", "can")

    #step 1: generate the data for the graph
    graph = {}
    for atom in ob.OBMolAtomIter(mol):
        a = atom.GetId()
        bonds = []
        for bond in ob.OBAtomBondIter(atom):
            b1 = bond.GetBeginAtom().GetId()
            b2 = bond.GetEndAtom().GetId()
            if b1 != a:
                bonds.append(b1)
            else:
                bonds.append(b2)
        graph[a] = bonds

    #step 2: generate all subgraphs of max_size
    subgraphs = []
    for atom in ob.OBMolAtomIter(mol):
        for subgraph in traverse(atom.GetId(),graph,max_size):
            if subgraph not in subgraphs:
                subgraphs.append(subgraph)

    #step 3: generate smiles for all subgraphs by deleting all other atoms
    for subgraph in subgraphs:
        #copy molecule
        copy = ob.OBMol(mol)
        #delete all atoms not in the subgraph
        for atom in ob.OBMolAtomIter(mol): #iterate over atoms of mol, NOT copy, because atoms are deleted in place
            idx = atom.GetId()
            if idx not in subgraph: #all atoms in subgraph are kept
                #increment implicit H counts of bonding partners
                for bond in ob.OBAtomBondIter(atom):
                    a1 = bond.GetBeginAtom().GetId()
                    a2 = bond.GetEndAtom().GetId()
                    bo = bond.GetBondOrder()
                    if a1 != idx:
                        if a1 in subgraph: #only kept atoms need to be handled
                            copy.GetAtomById(a1).SetImplicitHCount(copy.GetAtomById(a1).GetImplicitHCount()+bo)
                    else:
                        if a2 in subgraph: #only kept atoms need to be handled
                            copy.GetAtomById(a2).SetImplicitHCount(copy.GetAtomById(a2).GetImplicitHCount()+bo)
                copy.DeleteAtom(copy.GetAtomById(idx))

        #convert to smiles
        smiles = obConversion.WriteString(copy).strip()
        obConversion.ReadString(copy, smiles) #this is done to really get canonical smiles
        smiles = obConversion.WriteString(copy).strip()

        if smiles not in already_generated:
            already_generated.add(smiles)
            yield copy

