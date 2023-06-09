#Imports required packages
import os
import xml.etree.ElementTree as ET
import sys
import csv
import mygene 
import requests
from functools import cached_property, lru_cache
import urllib
import pandas as pd
import numpy as np
import traceback
from io import StringIO
from chembl_webresource_client.new_client import new_client
import re
mg = mygene.MyGeneInfo()

#Utility dictionaries containing element/mass mappings as well as FASTA mapping
FASTAdict = {'ALA':'A', 'ALX':'B', 'CYS': 'C', 'ASP' : 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS':'K', 'LEU':'L', 'MET': 'M', 'ASN':'N', 'PRO':'P', 'GLN': 'Q', 'ARG':'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'UNK':'X', 'TYR':'Y', 'GLX':'Z'}
atom_dict = {'H':1.01, 'C':12.01, 'O':16.00, 'N':14.01, 'P':30.97, 'F':18.998, 'S':32.06, 'B':10.81, 'K':39.1, 'I':126.904, 'BR':79.904, 'CL':35.453, 'CA':40.08, 'NA':22.99, 'MG':24.305, 'AL':26.98, 'CR':51.996, 'NE':20.179, 'BE':9.01, 'FE':55.847, 'CO':58.933,'AG':107.868, 'CD':112.41, 'NI':58.693, 'ZN':65.39, 'BE':9.0122, 'IN':114.818, 'SI':28.085, 'SC':44.956, 'TI':47.867, 'V':50.941, 'MN':54.938, 'CU':63.546, 'GA':59.723, 'GE':72.64, 'SE':78.96, 'KR':83.8, 'ZR':91.224, 'NB':92.906, 'PD':106.42, 'SN':118.71, 'SB':121.76, 'XE':131.293, 'BA':137.327, 'LA':138.91, 'LI':6.941, 'HG':200.59, 'PB':207.2, 'BI':208.98, 'PO':209, 'TI':204.3833, 'AU':196.9665, 'IR':192.217, 'PT':195.078, 'RE':186.207, 'W':183.84, 'TA':180.948, 'YB':173.04, 'EU':151.964, 'ND':144.25, 'CE':140.116, 'TH':232.04, 'U':238.029, 'PU':244, 'FR':223, 'PA':231.04, 'HO':164.93, 'SM':150.36, 'PR':140.908, 'TE':127.6, 'TC':98, 'Y':88.906}
atom_keys = atom_dict.keys()

class atom:
    """
    Atom class

    This class contains the essential parameters of an individual atom in a protein structural file, including the coordinates, element, atomic mass, parent residue, and line data.
    """
    def __init__(self, element, x, y, z, data, residue = ''):
        
        self.element = element
        self.x = x
        self.y = y
        self.z = z
        self.residue = residue
        self.data = data
        self.mass = atom_dict[self.element]
             
class residue:
    """
    Residue class
    
    This class contains the essential parameters of an individual residue in a protein structural file, including the child atoms, protein sidechain, residue index, and type of amino acid.
    """
    def __init__(self, amino_acid, index, chain):
        
        self.amino_acid = amino_acid
        self.index = index
        self.atoms = []
        self.chain = chain
            
    @cached_property
    def center(self):
        """
        Residue Center of Mass
    
        This property returns the center of mass of the residue, which is generated from all atom components. 
    
        Returns
        -------
        list
            A list consisting of the x, y, and z coordinate in 3D space of the residue center of mass position for the protein structural file.
        """
        x, y, z, totmass = 0, 0, 0, 0
        for atom in self.atoms:
            x += atom.x
            y += atom.y 
            z += atom.z
            totmass += atom.mass
        x /= atom.mass
        y /= atom.mass
        z /= atom.mass
        return [x, y, z]
        
    def __getitem__(self, key):
        """
        Residue Index
    
        This property returns the atom contained in the parent residue based on the index of the atom as it appears in the structural file.
        
        Parameters
        ----------
        index : int
            The atom index

        Returns
        -------
        class
            Returns the atom class contained at the given index.
        """
        try:  
            if isinstance(key, slice):
                atoms = [i for index, i in enumerate(self.atoms) if index > key.start and index <= key.stop] 
                return atoms
            else:
                return self.residue_list[key]
        except IndexError:
            print('Atom index out of bounds!')
        except:
            pass
            
    @lru_cache
    def __len__(self):
        """
        Residue Length
    
        This property returns the total number of atoms contained with the residue.
    
        Returns
        -------
        int
            A number of all the atoms contained within the specific residue.
        """
        return len(self.atoms)
        
class ligand:
    """
        Ligand Class
    
        This property returns the center of mass of the residue, which is generated from all atom components. 
        """
    def __init__(self, ID):
        self.ID = ID
        self.sites = []
        self.atoms = []
        
   
    def download(self, directory = os.getcwd()):
        """
        Download Ligand
    
        This method downloads the specific ligand structural file based on the ligand ID. 
    
        Parameters
        ----------
        arg1 : str
            Directory structural file will be downloaded in. Defaults to current user directory.
        """
        url = 'https://files.rcsb.org/ligands/download/' + self.ID + '_ideal.sdf'
        urllib.request.urlretrieve(url, directory + '/' + self.ID + '.sdf')
        
    @cached_property
    def center(self):
        """
        Ligand Center of Mass
    
        This property returns the center of mass of the ligand, which is generated from all atom components. 
    
        Returns
        -------
        list
            A list consisting of the x, y, and z coordinate in 3D space of the ligand center of mass position for the specific ligand pose registered within the file.
        """
        x, y, z, totmass = 0, 0, 0, 0
        for atom in self.atoms:
            x += atom.x
            y += atom.y
            z += atom.z
            totmass += atom.mass
        x /= totmass
        y /= totmass
        z /= totmass
        return [x, y, z]
        
#Main protein class
class Protein:
    """
    Protein Class

    This class represents an individual protein, it's provided ID and relevant structural file (generated via Alphafold or scraped from RCSB), and analytical methods and properties for biostatistics and bioinformatics purposes.
    """
    def __init__(self, protein_name, protein_type = 'PDB', species = 'human'):
        self.protein = protein_name
        self.ID_type = protein_type
        self.species = species
        
        #Accepts Uniprots if they are the correct length and user forgot to specify
        if self.ID_type == 'PDB' and len(self.protein) > 4:
            self.ID_type = 'Uniprot'
        
        #length of each side chain
        self.chains = {}
        
    def __len__(self):
        """
        Protein Length
    
        This property returns the total number of residues within the protein structural file.
        
        Returns
        -------
        int
            The total number of residues within the protein.
        """
        length = 0
        try:
            for chain in self.chains:
                length += self.chains[chain]
        except:
            length = 0
            print('Protein structural file not downloaded!')
            
        return length
    
    def __getitem__(self, key):
        """
        Residue Index
    
        This magic method allows the user to index the protein for individual residues by either slicing or accessing an individual residue.
    
        Parameters
        ----------
        slice : slice
            A slice of the residues within the protein file.
            
        key : int
            An index of the residue within the protein file.
    
        Returns
        -------
        class, list
        Returns the class or the list of classes queried by the index.
        """
        try:  
            if isinstance(key, slice):
                residues = [i for index, i in enumerate(self.residue_list) if index > key.start and index <= key.stop] 
                return residues
            else:
                return self.residue_list[key]
        except IndexError:
            print('Residue index out of bounds!')
        except:
            print(traceback.format_exc())
            print('Protein file has not been downloaded yet! Download the protein with <protname>.download()')
            )
    def download(self, destination = os.getcwd()):
        """
        Download Protein
    
        This method allows the user to download the simulated structural file of the specific protein from Alphafold (if a non-PDB ID is given) or a true structural file from RCSB (if a specific PDB ID is provided). Can identify process cif/mmcif files as well as .pdb files. Downloading the protein also collects all residue and atom information from the structural file, populating the attributes of the protein class.
    
        Parameters
        ----------
        destination : str
            The destination of the protein file, defaulting to the current user directory.
        """
        self.cif = False
                
        #Downloads files from RCSB if they are PDB files
        if self.ID_type == 'PDB':
            try:
                protname = self.protein + '.pdb'
                urllib.request.urlretrieve('http://files.rcsb.org/download/' + protname, destination + '/' + protname)
            except:
                try:
                    protname = self.protein + '.cif'
                    urllib.request.urlretrieve('http://files.rcsb.org/download/' + protname, destination + '/' + protname)
                    self.cif = True
                except:
                    print(traceback.format_exc())
                    print('Unable to identify valid PDB or CIF file for ' + self.protein + '!')
        elif self.ID_type.upper() == 'UNIPROT':
            try:
                url = 'https://alphafold.ebi.ac.uk/files/'
                for i in range(0, 9):
                    try:
                        alphaname = 'AF-' + self.protein + '-F1-model_v' + str(i) + '.pdb'
                        urllib.request.urlretrieve(url + alphaname, destination + '/' + alphaname) 
                        protname = self.protein + '.pdb'
                    except:
                        continue
                    os.rename(alphaname, protname)
            except:
                try:
                    for i in range(0, 9):
                        try:
                            alphaname = 'AF-' + self.protein + '-F1-model_v' + str(i) + '.pdb'
                            urllib.request.urlretrieve(url + alphaname, destination + '/' + alphaname) 
                            protname = self.protein + '.pdb'
                        except:
                            continue
                        os.rename(alphaname, protname)
                    self.cif = True
                except:
                    print('Unable to identify Alphafold structure for protein ' + self.protein + '!')
        elif self.ID_type.upper() == 'HGNC':
            try:
                out = mg.querymany(self.protein, scopes = 'symbol', fields = 'uniprot', species = self.species)
            except:
                out = mg.querymany(self.protein, scopes = 'symbol', fields = 'uniprot', species = 'all')
            try:
                uniprot = out[0]['uniprot']['Swiss-Prot']
                for i in range(0, 9):
                    try:
                        alphaname = 'AF-' + uniprot + '-F1-model_v' + str(i) + '.pdb'
                        urllib.request.urlretrieve(url + alphaname, destination + '/' + alphaname) 
                        protname = self.protein + '.pdb'
                    except:
                        continue
                    os.rename(alphaname, protname)
            except:
                try:
                    line = 'https://rest.uniprot.org/uniprotkb/search?query=' + self.protein + '%20' + self.species
                    Data = requests.get(line).text
                    firstResultIndex = Data.index('primaryAccession')
                    uniprot = Data[firstResultIndex + 19:firstResultIndex + 25]
                    for i in range(0, 9):
                        try:
                            alphaname = 'AF-' + uniprot + '-F1-model_v' + str(i) + '.pdb'
                            urllib.request.urlretrieve(url + alphaname, destination + '/' + alphaname) 
                            protname = self.protein + '.pdb'
                        except:
                            continue
                    os.rename(alphaname, protname)
                except:
                    print(self.protein + ' not found in any libraries!')
        else:
            try:
                line = 'https://rest.uniprot.org/uniprotkb/search?query=' + self.protein + '%20' + self.species
                Data = requests.get(line).text
                firstResultIndex = Data.index('primaryAccession')
                uniprot = Data[firstResultIndex + 19:firstResultIndex + 25]
                for i in range(0, 9):
                    try:
                        alphaname = 'AF-' + uniprot + '-F1-model_v' + str(i) + '.pdb'
                        urllib.request.urlretrieve(url + alphaname, destination + '/' + alphaname) 
                        protname = self.protein + '.pdb'
                    except:
                        continue
                    os.rename(alphaname, protname)
            except:
                try:
                    url = 'https://alphafold.ebi.ac.uk/files/'
                    for i in range(0, 9):
                        try:
                            alphaname = 'AF-' + uniprot + '-F1-model_v' + str(i) + '.pdb'
                            urllib.request.urlretrieve(url + alphaname, destination + '/' + alphaname) 
                            protname = self.protein + '.pdb'
                        except:
                            continue
                    os.rename(alphaname, protname)
                    self.cif = True
                except:
                    print(self.protein + ' not found in any libraries!')
        
        #Atomizes protein file
        if self.cif == False:
            with open(self.protein + '.pdb', 'r') as f:
                data = f.readlines()
            self.residue_list, curr_residue, curr_lig, currnum, curr_chain = [], -1, '', 0, ''
            if self.ID_type == 'PDB':
                current_ligand_index, self.ligand_list = -1, []
            for count, line in enumerate(data):
                if line[0:6] == 'SEQRES':
                    chain = line[11:12]
                    if curr_chain != chain:
                        curr_chain = chain
                        self.chains[chain] = int(line[13:17].strip())
                        
                elif line[0:4] == 'ATOM':
                    resnum = int(line[23:27].strip())
                    if resnum != curr_residue:
                        curr_residue = resnum
                        chain, AA = line[21:22], line[17:20]
                        res = residue(amino_acid = AA, index = resnum, chain = chain)
                        self.residue_list.append(res)
                    element = line[77:79].strip()
                    if element not in atom_keys:
                        element = element[:1]
                    atm = atom(element, float(line[30:38].strip()), float(line[38:47].strip()), float(line[47:56].strip()), residue = res.chain + str(resnum), data = line)
                    res.atoms.append(atm)
            
            #Generates ligands if PDB filetype and ligands in file
            if self.ID_type == 'PDB':
                for count, line in enumerate(data):
                    if line[0:6] == 'HETATM':
                        ligands = line[17:20]
                        if ligands == 'HOH':
                            break
                        if ligands != curr_lig:
                            lignum = int(line[22:27])
                            if lignum != currnum:
                                currnum = lignum
                                lig = ligand(ligands)
                                self.ligand_list.append(lig)
                                current_ligand_index += 1
                            element = line[77:79].strip()
                        if element not in atom_keys:
                            element = element[:1]
                        atm = atom(element, float(line[30:38].strip()), float(line[38:47].strip()), float(line[47:56].strip()), data = line)
                        self.ligand_list[current_ligand_index].atoms.append(atm)
        else:
            with open(self.protein + '.cif', 'r') as f:
                data = f.readlines()
            
            #Obtains chains from protein file
            numbers, letters, curr_num = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'], ['A', 'B',' C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'], 1
            for count, line in enumerate(data):
                if len(line) < 16:
                    if line[:1] in numbers:
                        line = [i for i in line.split(' ') if i != '']
                        chain_num = int(line[0])
                        if chain_num != curr_num:
                            line = [i for i in data[count - 1].split(' ') if i != '']
                            curr_num = chain_num
                            self.chains[letters[chain_num - 2]] = int(line[1])
                if len(self.chains) > 0 and line[:1] == '#':
                    chain_num = int(data[count - 1][:1])
                    line = [i for i in data[count - 1].split(' ') if i != '']
                    self.chains[letters[chain_num]] = int(line[1])
                    break
                    
            #Obtains residues from protein file
            self.residue_list, curr_residue = [], -1
            for count, line in enumerate(data):
                if line[0:4] == 'ATOM':
                    resnum = int(line[33:37].strip())
                    if resnum != curr_residue:
                        curr_residue = resnum
                        chain, AA = line[28:29], line[24:27]
                        res = residue(amino_acid = AA, index = resnum, chain = chain)
                        self.residue_list.append(res)
                    element = line[14:16].strip()
                    if element not in atom_keys:
                        element = element[:1]
                    atm = atom(element, float(line[39:47].strip()), float(line[47:55].strip()), float(line[55:63].strip()), residue = res.chain + str(resnum), data = line)
                    res.atoms.append(atm)
                    
        #Strips binding sites if available
        counter, firstRun, curr_site, sites = 0, True, '', []
        if self.ID_type == 'PDB' and self.cif == False:
            ligands = []
            for myligand in self.ligand_list:
                ligands.append(myligand.ID)
            for count, line in enumerate(data):
                if line[0:4] == 'SITE':
                    ID = line[11:14]
                    res = line.split(' ')
                    for number, site in enumerate(res):
                        site = line[22 + 11 * number:28 + 11 * number].strip()
                        if site != '' and '\n' not in site and line[18:21] != 'HOH':
                            sites.append(site)
                    if firstRun == True:
                        curr_site = line[11:14]
                        firstRun = False
                    if line[11:14] != curr_site or data[count + 1][0:4] != 'SITE':
                        curr_site = line[11:14]
                        self.ligand_list[counter].sites = sites
                        sites = []
                        counter += 1    
    
        #Generates FASTA sequence based on protein file
        self.FASTA = ''
        for id_residue in self.residue_list:
            self.FASTA += FASTAdict[id_residue.amino_acid]
    
    def to_csv(self, destination = os.getcwd()):
        """
        Atoms to CSV
    
        This method returns a .csv file of each individual atom of the protein and all contained attributes within the atom.
    
        Parameters
        ----------
        destination : str
            Destiation of the csv file. Defaults to the user's current directory
        """
        atoms_x, atoms_y, atoms_z, atoms_element, atoms_residue, residue_index, residue_chain = [], [], [], [], [], [], []
        for count, residue in enumerate(self.residue_list):
            for atom in residue.atoms:
                atoms_x.append(atom.x)
                atoms_y.append(atom.y)
                atoms_z.append(atom.z)
                atoms_element.append(atom.element)
                atoms_residue.append(residue.amino_acid)
                residue_index.append(residue.index)
                residue_chain.append(residue.chain)
        atoms = pd.DataFrame.from_dict({'Atom x coordinate':atoms_x, 'Atom y coordinate':atoms_y, 'Atom z coordinate':atoms_z, 'Atom element':atoms_element, 'Atom residue':atoms_residue, 'residue index':residue_index, 'residue chain': residue_chain})
        atoms.to_csv(destination + '/Atoms' + self.protein + '.csv')
        
    @cached_property
    def center(self):
        """
        Protein Center of Mass
    
        This property returns the center of mass of the protein, which is generated from all atom components of each residue in the protein.
    
        Returns
        -------
        list
            A list consisting of the x, y, and z coordinate in 3D space of the protein center of mass position for the protein structural file.
        """
        try:
            x, y, z, totmass = 0, 0, 0, 0
            for res in self.residue_list:
                for atom in res.atoms:
                    x += atom.x
                    y += atom.y
                    z += atom.z
                    totmass += atom.mass
            x /= totmass
            y /= totmass
            z /= totmass
            return [x, y, z]
        except:
            print('Protein not downloaded yet or downloaded incorrectly!')
    
    @cached_property
    def sites(self):
        """
        Docked Ligand Center of Mass
    
        This property returns the center of mass for each ligand position of each unique ligand type within the protein file.
    
        Returns
        -------
        dict
            A dictionary containing the unique ligand identifier within the protein file (the key) and a list of the center of mass of each indivual ligand pose (the values)
        """
        if self.ID_type != 'PDB':
            return 'Not a PDB file with valid binding sites!'
        else:
            sites = {}
            for ligand in self.ligand_list:
                if len(ligand.sites) > 0:
                    try:
                        sites[ligand.ID].append(ligand.sites)
                    except:
                        sites[ligand.ID] = [ligand.sites]
                else:
                    sites = 'No valid binding sites listed in PDB file!'
                    break
        return sites
    
    @lru_cache
    def residues(self, residue_index):
        """
        Residue Indexing
    
        This method allows for indexing residues with a numeric or a side chain-numeric (i.e. A244) identifier.
    
        Parameters
        ----------
        index : int, str
            Numeric or string index of the residue.
            
        Returns
        -------
        class
            Returns the residue class at the given index.
        """
        try:
            try:
                index = int(re.sub('[^0-9]', '', residue_index))
                residue = self.residue_list[index]
                if residue.index != index:
                    diff = residue.index - index
                    residue = self.residue_list[index - diff]
            except:
                residue = self.residue_list[residue_index]
    
            residue_items = {'chain':residue.chain, 'AA':residue.amino_acid, 'index':residue.index, 'atoms':residue.atoms}
            res = pd.Series(residue_items, index = ['chain', 'AA', 'index', 'atoms'])
            return res
        except:
            return 'Index larger than number of residues!'
           
    @lru_cache
    def atoms(self, atom_index):
        """
        Atom Indexing
    
        This method allows for indexing atoms with a numeric identifier.
    
        Parameters
        ----------
        index : int
            Numeric index of the atom. 
            
        Returns
        -------
        class
            Returns the atom class at the given index.
        """
        ac = 0
        for res in self.residue_list:
            for atom in res.atoms:
                if ac == atom_index:
                    atom_items = {'element':atom.element, 'x':atom.x, 'y':atom.y, 'z':atom.z, 'residue':atom.residue, 'data':atom.data}
                    atm = pd.Series(atom_items, index = ['element', 'x', 'y', 'z', 'residue', 'data'])
                    return atm
                ac += 1
        return 'Index is larger than total number of atoms!'
    
    def concat(self, start, stop, destination = os.getcwd()):
        """
        Residue Concatenation
    
        This method allows for the concatenation of residues into a sub-structural file containing only the indexed residues.
    
        Parameters
        ----------
        start : int
            Starting numeric index of the first residue.
            
        stop: int
            Final numeric index of the last residue in the sub-structure.
            
        destination: str
            Location of the downloaded sub-structure, defaults to the user's current directory.
        """
        residues, atoms = [i for i in self.residue_list if i >= start and i < stop], []
        for res in residues:
            for atom in res.atoms:
                atoms.append(atom.data)
                
        with open(destination + '/' + self.protein + str(start) + 'to' + str(stop) + '.pdb', 'w') as f:
            for line in atoms:
                f.write(line)
    
    @cached_property
    def Uniprot(self):
        """
        Uniprot ID
    
        This properties allows for easy conversion of protein IDs to Uniprot IDs.

        Returns
        -------
        str
            Returns the Uniprot ID for the protein class.
        """
        if self.ID_type.upper() == 'UNIPROT':
            return self.protein
        elif self.ID_type.upper() == 'PDB':
            try:
                out = mg.querymany(self.protein, scopes = 'pdb', fields = 'uniprot', species = self.species)
            except:
                out = mg.querymany(self.protein, scopes = 'pdb', fields = 'uniprot', species = 'all')
            try:
                uniprot = out[0]['uniprot']
                try:
                    uniprot = out[0]['uniprot']['Swiss-Prot']
                except:
                    pass
                if uniprot != None:
                    return uniprot
                else:
                    raise KeyError
            except:
                try:
                    url = 'https://rest.uniprot.org/uniprotkb/search?query=' + self.protein + '%20' + self.species
                    Data = requests.get(url).text
                    firstResultIndex = Data.index('primaryAccession')
                    uniprot = Data[firstResultIndex + 19:firstResultIndex + 25] 
                
                    return uniprot
                except Exception as e:
                    print('Cannot find valid Uniprot ID for protein ' + self.protein + '!')
                    return 'N/A'
        elif self.ID_type.upper() == 'HGNC':
            try:
                out = mg.querymany(self.protein, scopes = 'symbol', fields = 'uniprot', species = self.species)
            except:
                out = mg.querymany(self.protein, scopes = 'symbol', fields = 'uniprot', species = 'all')
            try:
                uniprot = out[0]['uniprot']
                try:
                    uniprot = out[0]['uniprot']['Swiss-Prot']
                except:
                    pass
                if uniprot != None:
                    return uniprot
                else:
                    raise KeyError
            except:
               try:
                   url = 'https://rest.uniprot.org/uniprotkb/search?query=' + self.protein + '%20' + self.species
                   Data = requests.get(url).text
                   firstResultIndex = Data.index('primaryAccession')
                   uniprot = Data[firstResultIndex + 19:firstResultIndex + 25] 
                   return uniprot
               except:
                   print('Cannot find valid Uniprot ID for protein ' + self.protein + '!')
                   return 'N/A'
        else:
            try:
                url = 'https://rest.uniprot.org/uniprotkb/search?query=' + self.protein + '%20' + self.species
                Data = requests.get(url).text
                firstResultIndex = Data.index('primaryAccession')
                uniprot = Data[firstResultIndex + 19:firstResultIndex + 25] 
                return uniprot
            except:
                print('Cannot find valid Uniprot ID for protein ' + self.protein + '!')
                return 'N/A'
    
    @cached_property 
    def PDB(self):
        """
        PDB ID
    
        This properties allows for easy conversion of protein IDs to PDB IDs. A pandas dataframe containing the number of unique ligands, number of residues, and the resolution of the electron microscope used to image the protein is returned.

        Returns
        -------
        df
            Returns a dataframe containing relevant properties for all PDB IDs associated with the protein ID.
        """
        def scrape_data(PDB):  
            resdf = pd.DataFrame()
            resolutions, residue_number, lignums, proteins = [], [], [], []
            for ID in PDB:
                curr_lig, lignum = '', 0
                protname = ID + '.pdb'
                try:
                    urllib.request.urlretrieve('http://files.rcsb.org/download/' + protname, protname)
                except:
                    continue
                with open(protname, 'r') as f:
                    data, resolution = f.readlines(), 0
                for line in data:
                    if 'RESOLUTION' in line and line[:6] == 'REMARK':  
                        line = [l for l in line.split(' ') if l != '']
                        try:
                            resolution = float(line[3])
                        except:
                            break
                        break
                resolutions.append(resolution)
                curr_chain, tot_res = '', 0
                for count, line in enumerate(data):
                    if line[0:6] == 'SEQRES':
                        chain = line[11:12]
                        if curr_chain != chain:
                            curr_chain = chain
                            tot_res += int(line[13:17].strip())
                
                    if line[0:6] == 'HETATM':
                        ligands = line[17:20]
                        if ligands == 'HOH':
                            break
                        if ligands != curr_lig:
                            curr_lig = ligands
                            lignum += 1
                lignums.append(lignum)    
                residue_number.append(tot_res)
                os.remove(protname)
                proteins.append(ID)
                
            df = pd.DataFrame.from_dict({'PDB': proteins, 'Resolution':resolutions, 'Residue Number':residue_number, 'Unique Ligands': lignums})
            df = df[df['Resolution'] > 0]
            df = df.sort_values(by = 'Resolution', ascending = True, ignore_index = True)
            return df         
        if self.ID_type == 'PDB': 
            return scrape_data([self.protein])
        elif self.ID_type == 'Uniprot':
            try:
                try:
                    out = mg.querymany(self.protein, scopes = 'uniprot', fields = 'pdb', species = self.species)
                except:
                    out = mg.querymany(self.protein, scopes = 'uniprot', fields = 'pdb', species = 'all')
                pdb = out[0]['pdb']
                return scrape_data(pdb)
            except:
                try:
                    PDB_IDs, index = [], 0
                    try:
                        url = 'https://rest.uniprot.org/uniprotkb/search?query=' + self.protein + '%20' + self.species
                        Data = requests.get(url).text
                        while True:
                            firstResultIndex = Data.index('"PDB"')
                            pdb = Data[firstResultIndex + 12:firstResultIndex + 16]
                            Data =  Data[firstResultIndex + 20:]
                            PDB_IDs.append(pdb)
                    except:
                        pass
                    PDB_IDs = list(set(PDB_IDs))
                    return scrape_data(PDB_IDs)
                except:
                    print('Cannot find valid PDB ID for protein ' + self.protein + '!')
                    return 'N/A'
        elif self.ID_type == 'HGNC':
            try:
                try:
                    out = mg.querymany(self.protein, scopes = 'symbol', fields = 'pdb', species = self.species)
                except:
                    out = mg.querymany(self.protein, scopes = 'symbol', fields = 'pdb', species = 'all')
                pdb = out[0]['pdb']
                return scrape_data(pdb)
            except:
                try:
                    PDB_IDs, index = [], 0
                    try:
                        url = 'https://rest.uniprot.org/uniprotkb/search?query=' + self.protein + '%20' + self.species
                        Data = requests.get(url).text
                        while True:
                            firstResultIndex = Data.index('"PDB"')
                            pdb = Data[firstResultIndex + 12:firstResultIndex + 16]
                            Data =  Data[firstResultIndex + 20:]
                            PDB_IDs.append(pdb)
                    except:
                        pass
                    return scrape_data(PDB_IDs)
                except:
                    print('Cannot find valid PDB ID for protein ' + self.protein + '!')
                    return 'N/A'
        else:
            try:
                PDB_IDs, index = [], 0
                try:
                    url = 'https://rest.uniprot.org/uniprotkb/search?query=' + self.protein + '%20' + self.species
                    Data = requests.get(url).text
                    while True:
                        firstResultIndex = Data.index('"PDB"')
                        pdb = Data[firstResultIndex + 12:firstResultIndex + 16]
                        Data =  Data[firstResultIndex + 20:]
                        PDB_IDs.append(pdb)
                except:
                    pass
                return scrape_data(PDB_IDs)
            except:
                print('Cannot find valid PDB ID for protein ' + self.protein + '!')
                return 'N/A'
        
    
    @cached_property
    def Gene(self):
        """
        Gene ID
    
        This properties allows for easy conversion of protein IDs to Gene IDs.

        Returns
        -------
        str
            Returns the Gene ID for the protein class.
        """
        def search_uniprot(protein):
            url = 'https://rest.uniprot.org/uniprotkb/search?query=' + self.protein + '%20' + self.species
            Data = requests.get(url).text
            firstResultIndex = Data.index('geneName":')
            Gene = Data[firstResultIndex + 20:firstResultIndex + 31]
            Gene = Gene.split('"')[0]
            return Gene
            
        if self.ID_type == 'Gene':
            return self.protein
        elif self.ID_type == 'Uniprot':
            try:
                try:
                    out = mg.querymany(self.protein, scopes = 'uniprot', fields = 'symbol', species = self.species)
                except:
                    out = mg.querymany(self.protein, scopes = 'uniprot', fields = 'symbol', species = 'all')
                gene = out[0]['symbol']
                return gene
            except:
                return search_uniprot(self.protein)
        elif self.ID_type == 'PDB':
             try:
                 try:
                     out = mg.querymany(self.protein, scopes = 'pdb', fields = 'symbol', species = self.species)
                 except:
                     out = mg.querymany(self.protein, scopes = 'pdb', fields = 'symbol', species = 'all')
                 gene = out[0]['symbol']
                 return gene
             except:
                return search_uniprot(self.protein)
            
    @cached_property
    def ChEMBL(self):
        """
        ChEMBL ID
    
        This properties allows for easy conversion of protein IDs to ChEMBL IDs.

        Returns
        -------
        str
            Returns the ChEMBL ID for the protein class.
        """
        try:
            target = new_client.target
            gene_name = self.Gene
            res = target.filter(target_synonym__icontains=gene_name).only('target_chembl_id')[0]
            ChEMBL = res['target_chembl_id']
            return ChEMBL
        except:
            try:
                url = 'https://rest.uniprot.org/uniprotkb/search?query=' + self.protein + '%20' + self.species
                Data = requests.get(url).text
                firstResultIndex = Data.index('ChEMBL')
                ChEMBL = Data[firstResultIndex + 14:firstResultIndex + 34]
                ChEMBL = ChEMBL.split('"')
                return ChEMBL[0]
            except:
                print('No ChEMBL ID found for ' + self.protein + '!')
                return 'N/A'
         
    @cached_property
    def ligands(self):
        """
        Protein Ligands
    
        This property returns a series containing the ID of the primary ligand contained within the file, as well as the ID of all present confactors. Only returns ligands if the file is a PDB file with valid ligand structures contained.
    
        Returns
        -------
        df
            A pandas dataframe containing the primary ligand IDs and cofactor IDs within the protein structural file.
        """
        ligs = []
        for ligand in self.ligand_list:
            if ligand.ID not in ligs:
                ligs.append(ligand.ID)
            ligand = ligs[1]
            cofactors = ligs[1:]
        ligands = pd.DataFrame.from_dict({'Primary Ligand':ligand, 'Cofactors':cofactors})
        return ligands
    
    @cached_property   
    def interactions(self):
        """
        Protein Interactions
    
        This property allows the user to query the Bioinformatics databases STITCH, BindingDB, and ChEMBL for relevant protein-ligand and protein-protein interactions. BindingDB interactions are contained in a dataframe containing the Kd constant in nanoMolar. ChEMBL IDs are returned along with the IC50 or DC50 values in nanoMolar. STITCH data represents a dictionary-embedded list of all proteins and ligands with network-based interactions. If valid data cannot be pulled from any site, their output datatypes will not be included in output.
    
        Returns
        -------
        list
            A list containing two dataframes for BindingDB and ChEMBL protein-ligand intearction data, as well as a dictionary containing all STITCH interactors for the given protein ID with a high confidence of interaction.
        """
        try:
            #Queries bindingDB for ligand interactions
            try:
                Affinities, ligands = [], []
                url = 'https://bindingdb.org/axis2/services/BDBService/getLigandsByUniprot?uniprot=' + self.Uniprot
                mytree = ET.fromstring(requests.get(url, timeout=4).text)
                for tree in mytree:
                    if 'affinities' in i.tree:
                        ligands.append(i[1].text)
                        Affinities.append(float(i[3].text) * 1000)
                BindingDB = pd.DataFrame.from_dict({'Ligands':ligands, 'Kd': Affinities})
            except:
                print('BindingDB request timed out!')
                
            #Queries STITCH for protein and ligand interactions
            try:
                url = 'http://string-db.org/api/tsv-no-header/resolve?identifier=' + self.Uniprot
                Result = [i for i in re.split(' |\t', requests.get(url).text.strip()) if i != '']
                if Result[3].upper() == Result[3]:
                    Gene = Result[3]
                else:
                    Gene = Result[4]
                ID = Result[0]
                url = 'http://string-db.org/api/tsv-no-header/interactorsList?identifiers=' + ID + '&required_score=900'
                tsv = pd.read_table(StringIO(requests.get(url).text), names = ['Conec1', 'Conec2', 'ID', 'Gene', 'CID', 'Prob', 'Z', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6'], header = None)
                IDlist = []
                for count, row in tsv.iterrows():
                    ID = row['Gene']
                    if ID != Gene:
                        IDlist.append(ID)
                STITCH = {'STITCH':set(IDlist)}
            except:
                print('No compounds or proteins found in STICH!')
                
            #Queries ChEMBL for ligand interactions
            try:
                target, ligands, activities = new_client.target, [], []
                activity = new_client.activity
                activities = activity.filter(target_chembl_id=self.ChEMBL).filter(standard_type="IC50")
                for act in activities:
                    value = float(act['value'])
                    if act['units'] != 'nM':
                        if act['units'] == 'uM':
                            value *= 1000
                        elif act['units'] == 'pM':
                            value /= 1000
                    try:
                        ligands.append(act['canonical_smiles'])
                        activities.append(value)
                    except:
                        continue
                ChEMBL = pd.DataFrame.from_dict({'Ligands':ligands, 'IC50 (uM)':activities})
                ChEMBL = {'ChEMBL':chembl}
            except:
                print('No compounds or proteins found in ChEMBL!')
            results = []  
            try:
                results.append(BindingDB)
            except:
                pass
            try:
                results.append(ChEMBL)
            except:
                pass
            try:
                results.append(STITCH)
            except:
                pass
            return results
        except:
            print(traceback.format_exc())
            return 'N/A'
    '''
    TODO
    - Add secondary structure prediction to residues (use s4pred?)
    - Add more ID conversions
    '''
'''
#For unit testing    
if __name__ == '__main__':
    p = Protein('1HK4')
    p.download()
    print(p.ligands)
    #new = Protein(p.Uniprot)
    #new.download()
    #print(new.interactions)
    #print(protein.interactions)
'''