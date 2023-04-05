#Imports required packages
from import_modules import *
mg = mygene.MyGeneInfo()

#Utility files containing element/mass mappings as well as FASTA mapping
from tiger_utils import *

#Subclass for individual atoms
class atom:
    def __init__(self, element, x, y, z, data, residue = ''):
        self.element = element
        self.x = x
        self.y = y
        self.z = z
        self.residue = residue
        self.data = data
        self.mass = atom_dict[self.element]
             
#Subclass for individual protein residues
class residue:
    def __init__(self, amino_acid, index, chain):
        self.amino_acid = amino_acid
        self.index = index
        self.atoms = []
        self.chain = chain
            
    #Calculates center of residue
    def calc_center(self):
        x, y, z, totmass = 0, 0, 0, 0
        for atom in self.atoms:
            x += atom.x
            y += atom.y 
            z += atom.z
        
    #Allows user to access residue atoms via indexing
    def __getitem__(self, key):
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
            
    #Returns number of atoms in residue
    @lru_cache
    def __len__(self):
        return len(self.atoms)
        
#Subclass for ligands in protein file (if present)
class ligand:
    def __init__(self, ID):
        self.ID = ID
        self.sites = []
        self.atoms = []
        
    #Downloads ligand from PDB database
    def download_ligand(self, directory = os.getcwd()):
        url = 'https://files.rcsb.org/ligands/download/' + self.ID + '_ideal.sdf'
        urllib.request.urlretrieve(url, directory + '/' + self.ID + '.sdf')
            
#Main protein class
class Protein:
    def __init__(self, protein_name, protein_type = 'PDB', species = 'human'):
        self.protein = protein_name
        self.ID_type = protein_type
        self.species = species
            
        #length of each side chain
        self.chains = {}
        
    #Returns length of protein
    def __len__(self):
        length = 0
        for chain in self.chains:
            length += self.chains[chain]
        return length
    
    #Allows user to access residues via indexing
    def __getitem__(self, key):
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
            
    #Attempts to download protein file to a default directory (or a user-specified one)
    def download(self, destination = os.getcwd()):
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
                alphaname = 'AF-' + self.protein + '-F1-model_v2.pdb'
                request.urlretrieve(url + alphaname, destination + '/' + alphaname) 
                protname = self.protein + '.pdb'
                os.rename(alphaname, protname)
            except:
                try:
                    url = 'https://alphafold.ebi.ac.uk/files/'
                    alphaname = 'AF-' + self.protein + '-F1-model_v2.cif'
                    request.urlretrieve(url + alphaname, destination + '/' + alphaname) 
                    protname = self.protein + '.cif'
                    os.rename(alphaname, protname)
                    self.cif = True
                except:
                    print('Unable to identify Alphafold structure for protein ' + self.protein + '!')
        elif self.ID_type.upper() == 'HGNC':
            out = mg.querymany(self.protein, scopes = 'symbol', fields = 'uniprot', species = self.species)
            try:
                uniprot = out[0]['uniprot']['Swiss-Prot']
                url = 'https://alphafold.ebi.ac.uk/files/'
                alphaname = 'AF-' + uniprot + '-F1-model_v2.pdb'
                request.urlretrieve(url + alphaname, destination + '/' + alphaname)  
                protname = self.protein + '.pdb'
                os.rename(alphaname, protname)
            except:
                try:
                    line = 'https://rest.uniprot.org/uniprotkb/search?query=' + self.protein + '%20' + self.species
                    Data = requests.get(line).text
                    firstResultIndex = Data.index('primaryAccession')
                    uniprot = Data[firstResultIndex + 19:firstResultIndex + 25]
                    url = 'https://alphafold.ebi.ac.uk/files/'
                    alphaname = 'AF-' + uniprot + '-F1-model_v2.pdb'
                    request.urlretrieve(url + alphaname, destination + '/' + alphaname)  
                    protname = self.protein + '.pdb'
                    os.rename(alphaname, protname)
                except:
                    print(self.protein + ' not found in any libraries!')
        else:
            try:
                line = 'https://rest.uniprot.org/uniprotkb/search?query=' + self.protein + '%20' + self.species
                Data = requests.get(line).text
                firstResultIndex = Data.index('primaryAccession')
                uniprot = Data[firstResultIndex + 19:firstResultIndex + 25]
                url = 'https://alphafold.ebi.ac.uk/files/'
                alphaname = 'AF-' + uniprot + '-F1-model_v2.pdb'
                request.urlretrieve(url + alphaname, destination + '/' + alphaname)  
                protname = self.protein + '.pdb'
                os.rename(alphaname, protname)
            except:
                try:
                    url = 'https://alphafold.ebi.ac.uk/files/'
                    alphaname = 'AF-' + uniprot + '-F1-model_v2.cif'
                    request.urlretrieve(url + alphaname, destination + '/' + alphaname) 
                    protname = self.protein + '.cif'
                    os.rename(alphaname, protname)
                    self.cif = True
                except:
                    print(self.protein + ' not found in any libraries!')
        
        #Atomizes protein file
        if self.cif == False:
            with open(self.protein + '.pdb', 'r') as f:
                data = f.readlines()
            self.residue_list, curr_residue = [], -1
            if self.ID_type == 'PDB':
                curr_lig, currnum, curr_chain = '', 0, ''
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
    
    #Outputs all atoms in a .csv file
    def to_csv(self, destination = os.getcwd()):
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
        
    #Property for obtaining list of ligand sites
    @cached_property
    def sites(self):
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
    
    #Method for indexing protein residues
    @lru_cache
    def residues(self, residue_index):
        try:
            index = int(residue_index[1:])
            residue = self.residue_list[index]
            if residue.index != index:
                diff = residue.index - index
                residue = self.residue_list[index - diff]
        except:
            residue = self.residue_list[residue_index]

        residue_items = {'Chain':residue.chain, 'Amino Acid':residue.amino_acid, 'Index':residue.index, 'Atoms':residue.atoms}
        res = pd.Series(residue_items, index = ['Chain', 'Amino Acid', 'Index', 'Atoms'])
        return res
           
    #Method for indexing protein atoms
    @lru_cache
    def atoms(self, atom_index):
        ac = 0
        for res in self.residue_list:
            for atom in res.atoms:
                if ac == atom_index:
                    atom_items = {'Element':atom.element, 'x':atom.x, 'y':atom.y, 'z':atom.z, 'Residue':atom.residue, 'Line':atom.data}
                    atm = pd.Series(atom_items, index = ['Element', 'x', 'y', 'z', 'Residue', 'Line'])
                    return atm
                ac += 1
        return 'Index higher than total number of atoms!'
    
    #Method for reconstructing section of protein
    def concat(self, start, stop, destination = os.getcwd()):
        residues, atoms = [i for i in self.residue_list if i >= start and i < stop], []
        for res in residues:
            for atom in res.atoms:
                atoms.append(atom.data)
                
        with open(destination + '/' + self.protein + str(start) + 'to' + str(stop) + '.pdb', 'w') as f:
            for line in atoms:
                f.write(line)
    
    #Returns a uniprot identifier
    @cached_property
    def Uniprot(self):
        if self.ID_type.upper() == 'UNIPROT':
            return self.protein
        elif self.ID_type.upper() == 'PDB':
            out = mg.querymany(self.protein, scopes = 'pdb', fields = 'uniprot', species = self.species)
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
                   Data = requests.get(line).text
                   firstResultIndex = Data.index('primaryAccession')
                   uniprot = Data[firstResultIndex + 19:firstResultIndex + 25] 
                   return uniprot
               except:
                   print('Cannot find valid Uniprot ID for protein ' + self.protein + '!')
                   return 'N/A'
    
    #Returns a PDB identifier
    @cached_property 
    def PDB(self):
        if self.ID_type == 'PDB': 
            return self.protein
        elif self.ID_type.upper() == 'UNIPROT':
            out = mg.querymany(self.protein, scopes = 'uniprot', fields = 'pdb', species = self.species)
            try:
                pdb = out[0]['pdb']
                return pdb
            except:
                try:
                    url = 'https://rest.uniprot.org/uniprotkb/search?query=' + self.protein + '%20' + self.species
                    Data = requests.get(line).text
                    firstResultIndex = Data.index('PDB')
                    pdb = Data[firstResultIndex + 11:firstResultIndex + 15]
                    return pdb
                except:
                    print('Cannot find valid PDB ID for protein ' + self.protein + '!')
                    return 'N/A'
        elif self.ID_type.upper() == 'HGNC':
            out = mg.querymany(self.protein, scopes = 'uniprot', fields = 'pdb', species = self.species)
            try:
                pdb = out[0]['pdb']
                return pdb
            except:
                try:
                    url = 'https://rest.uniprot.org/uniprotkb/search?query=' + self.protein + '%20' + self.species
                    Data = requests.get(line).text
                    firstResultIndex = Data.index('PDB')
                    pdb = Data[firstResultIndex + 11:firstResultIndex + 15]
                    return pdb
                except:
                    print('Cannot find valid PDB ID for protein ' + self.protein + '!')
                    return 'N/A'
    
    #Returns HGNC identifier
    @cached_property
    def HGNC(self):
        if self.ID_type == 'HGNC':
            return self.protein
        else:
            url = 'https://rest.uniprot.org/uniprotkb/search?query=' + self.protein + '%20' + self.species
            Data = requests.get(url).text
            firstResultIndex = Data.index('HGNC:')
            HGNC = Data[firstResultIndex + 5:firstResultIndex + 15]
            HGNC = HGNC.split('"')[0]
            return HGNC
            
    #Returns ChEMBL identifier
    @cached_property
    def ChEMBL(self):
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
         
    #Property for obtaining protein binders. Searches ChEMBL, BindingDB, and STITCH
    @cached_property   
    def interactions(self):
        Affinities = {}
        try:
            #Queries bindingDB for ligand interactions
            try:
                url = 'https://bindingdb.org/axis2/services/BDBService/getLigandsByUniprot?uniprot=' + self.Uniprot
                mytree = ET.fromstring(requests.get(url, timeout=4).text)
                for tree in mytree:
                    if 'affinities' in i.tree:
                        Affinities[i[1].text] = float(i[3].text)
                BindingDB = {'BindingDB':Affinities}
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
                target, chembl = new_client.target, {}
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
                        temp = chembl[act['canonical_smiles']]
                        
                        chembl[act['canonical_smiles']].append(value)
                    except:
                        chembl[act['canonical_smiles']] = value
                ChEMBL = {'ChEMBL':chembl}
            except:
                print(traceback.format_exc())
                print('No compounds or proteins found in ChEMBL!')
             
            interaction_dict = {}
            try:
                interaction_dict.update(STITCH)
            except:
                pass
            try:
                interaction_dict.update(ChEMBL)
            except:
                pass
            try:
                interaction_dict.update(BindingDB)
            except:
                pass
            return interaction_dict
        except:
            print(traceback.format_exc())
            return 'N/A'
    '''
    TODO
    - Add secondary structure prediction to residues (use s4pred?)
    - Add more ID conversions
    - Add more identifirs from Uniprot (more cached properties?)
    
    FUTURE PROJECTS
    - Create documentation
        - Troubleshoot readthedocs (ask Abhi?)
    - Create ML software based around identifying similar sites on a certain residue in other proteins
        - Features are atom element, atom proximity, atom residue...
    '''
                    
#For unit testing    
if __name__ == '__main__':
    protein = Protein('1HK4')
    protein.download()
    print(protein.FASTA)
    #print(protein.residues('A44'))
    #print(len(protein))
    #print(protein[1:3])
    #print(protein.atoms(1000))
    #print(protein.interactions)