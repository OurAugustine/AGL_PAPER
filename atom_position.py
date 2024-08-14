#reading atom name and position
import re
import numpy as np

#extract postion from the pdb file
def extract_atom_positions_pdb(pdb_file: str ):
    atom_name_positions = []
    atom_cordinates = []
    with open(pdb_file, 'r') as f:
        intrest_atom=["C","N","O","S"]
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() in intrest_atom :
                atom_name = line[12:16].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                atom_name_positions.append((atom_name, (x, y, z)))
                atom_cordinates.append((x, y, z))
    return atom_name_positions ,atom_cordinates


 #extracting postion from the sdf file   
def extract_atom_positions_sdf(sdf_file : str)-> list:
    atom_name_positions = []
    atom_cordinates = []
    with open(sdf_file, 'r') as f:
        intrest_atom=["H","C","N","O","S","P","F","Cl","Br"]
        for line in f:
            if line.strip() == '$$$$':
                break  # End of molecule record
            elif len(line.strip().split()) == 9 and line.strip().split()[3] in intrest_atom:
                atom_name = line.strip().split()[3]
                x = float(line.strip().split()[0])
                y = float(line.strip().split()[1])
                z = float(line.strip().split()[2])
                atom_name_positions.append((atom_name, (x, y, z)))
                atom_cordinates.append((x, y, z))

    return atom_name_positions ,atom_cordinates

def siever_fxn(pro_atoms:list,lig_atoms:list,cut_off=9.0)->list:
    """this function takes list of atoms from protein and ligand the returs atoms that are only 9.0 A away from the ligand nearesr atom"""
    close_pro_atoms= []
    for pro_atom in pro_atoms:
        name,pro_cord = pro_atom
        pro_far=[]
        for lig_atom in lig_atoms:
            name,lig_cord=lig_atom
            dist = np.linalg.norm(np.array(pro_cord) - np.array(lig_cord))
            pro_far.append(dist)
        #min_dist=0
        if len(pro_far) !=0:
            min_dist= np.min(pro_far)
        else:
            min_dist=0

        if min_dist < cut_off:
            close_pro_atoms.append(pro_atom)
    return close_pro_atoms, lig_atoms


#extracting postion as cordinates , number of atoms and bonds from the sdf file   

def parse_sdf(file_path):
#this function will return atom cordinates ,count of attoms and bonds.It is the same as the second function
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    atoms = []
    bonds = []
    num_atoms = 0
    num_bonds = 0
    
    for i, line in enumerate(lines):
        if i == 3:  # This line contains the counts line
            counts_line = line.split()
            num_atoms = int(counts_line[0])
            num_bonds = int(counts_line[1])
        
        if 4 <= i < 4 + num_atoms:
            parts = line.split()
            x, y, z = map(float, parts[0:3])
            atoms.append((x, y, z))
        
        if 4 + num_atoms <= i < 4 + num_atoms + num_bonds:
            parts = line.split()
            atom1 = int(parts[0]) - 1  # Atom indices in SDF are 1-based
            atom2 = int(parts[1]) - 1
            bonds.append((atom1, atom2))
    #return the cordinates as atoms , bonds and number of atomsS
    return atoms, bonds, num_atoms

def main():
    #testing the first function
    pdb_file = "datasets/pdbbind_v2007/v2007/1a0q/1a0q_protein.pdb"
    atom_name_positions_pro,cords = extract_atom_positions_pdb(pdb_file)

    a=0
    for atom_name, position in atom_name_positions_pro:
        if a<10:#to print only first 10
            print(f'Atom: {atom_name}, Position: {position}')
        a=a+1
    count=0
    for position in cords:
        if count<10:#to print only first 10
            print(f'Position: {position}')
        count=count+1
    print('fucntion 1 success')
    #testing the second function
    sdf_file = "datasets/pdbbind_v2007/v2007/1a0q/1a0q_ligand.sdf"
    atom_name_positions_lig,cords = extract_atom_positions_sdf(sdf_file)
    count2=0
    for atom_name, position in atom_name_positions_lig:
        if count2<10:
            print(f'Atom: {atom_name}, Position: {position}')
        count2=count2+1
    print('function 2 succcess')
    print('')
    #testin the 3rd function which do the same work as function 2
    print('now checking if the parse function was a success')
    sdf_file = "datasets/pdbbind_v2007/v2007/1a0q/1a0q_ligand.sdf" 
    atoms, bonds, num_atoms = parse_sdf(sdf_file)
    print(num_atoms)
    close_pro,lig=siever_fxn(atom_name_positions_pro,atom_name_positions_lig,cut_off=9.5)
    print(len(atom_name_positions_pro))
    print(f'n.o of ligad atoms after : {len(lig)},original n.o of atoms : {len(atom_name_positions_lig)}')
    print(f'close proteins : {len(close_pro)},original n.o of atoms : {len(atom_name_positions_pro)}')    

if __name__ == "__main__":
    main()
