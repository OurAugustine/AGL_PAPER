import math
import numpy as np
from atom_position import extract_atom_positions_pdb # return name, cords
from atom_position import extract_atom_positions_sdf #returns name, cords
from atom_position import siever_fxn
from read_files import get_path_list
from ADJ_LAP_MAT  import create_combined_adjacency_matrix , create_laplacian_matrix ,fetch_eig_stat ,features_extractor
from Proteinfunctions import sep_core_non_dup,get_data_dict

#fetching the necessart protein name ,which are still subfolder names
ref_path="datasets/pdbbind_v2007/v2007/INDEX.2007.refined.data"  
core_path="datasets/pdbbind_v2007/v2007/INDEX.2007.core.data"
refined_folder_list,non_dup_folder_list,core_folder_list=sep_core_non_dup(ref_path,core_path)#returns alist of:all refined atoms,refined but non duplicated atoms

# Creating path to all protein function and ligand ,thes store them in alist of tuple 
main_folder = "datasets/pdbbind_v2007/v2007" 
extention1 = 'n.pdb'  # replace with your file extension
extention2= '.sdf'
paths_train = get_path_list(main_folder,refined_folder_list,'n.pdb','.sdf')#trainging sample paths
paths_test = get_path_list(main_folder,core_folder_list,'n.pdb','.sdf')
print("Number of expected items in the training path:",len(paths_train))
print("Number of expected items in the testin paths:",len(paths_test))

z=0
features_list=[] 
for item in paths_test:
    if z<5: #this if statemament is just to ensure we run only few of them as all will take alot of time
        
        pdb_file = item[0] 
        sdf_file = item[1]
        #print(pdb_file)
        #print(sdf_file)  

        pro_names_cord  , protein_atoms = extract_atom_positions_pdb (pdb_file)# getting cords and store in protein_atoms
        #ligand_atoms, _, _ = parse_sdf(sdf_file)#getting cords and store in ligad_atoms
        lig_names_cord,ligand_atoms=extract_atom_positions_sdf(sdf_file) #getting cords and store in ligad_atoms
        cutoff_distance = 9.1  # Defining a suitable cutoff distance
        

        #adjacency_matrix = create_combined_adjacency_matrix(protein_atoms, ligand_atoms, cutoff_distance)
        #laplacian_matrix = create_laplacian_matrix(adjacency_matrix)
        #features=fetch_eig_stat(laplacian_matrix)
        #features_list.append(features)
        #print(features)

        close_pro,lig = siever_fxn(pro_names_cord,lig_names_cord,cutoff_distance)
        features = features_extractor(close_pro,lig)
        features_list.append(features)
        
        print(len(features))
        print(len(features_list))


    z=z+1
print(features_list)

#remove the comment and run this lower part for complete extraction of X_data and Y_data for both training and test

""""
#creating y_train and y_test 
refined_proteins_dict=get_data_dict(ref_path)
core_proteins_dict=get_data_dict(core_path)
# List of keys
keys_ref = refined_folder_list
keys_core = core_folder_list

# Accessing values
y_tain = [float(refined_proteins_dict[key]) for key in keys_ref]
y_test = [float(core_proteins_dict[key]) for key in keys_core]

print(len(y_tain))  
print(len(y_test))


#creating X_train
X_train=[] 
for item in paths_train:
    pdb_file = item[0] 
    sdf_file = item[1]
    #print(pdb_file)
    #print(sdf_file)  

    pro_names_cord  , protein_atoms = extract_atom_positions_pdb (pdb_file)# getting cords and store in protein_atoms
    lig_names_cord,ligand_atoms=extract_atom_positions_sdf(sdf_file) #getting cords and store in ligad_atoms
    cutoff_distance = 9.1  # Defining a suitable cutoff distance

    close_pro,lig = siever_fxn(pro_names_cord,lig_names_cord,cutoff_distance)
    features = features_extractor(close_pro,lig)
    X_train.append(features)
        
    print(len(features))
    print(len(X_train))
#print(X_train)

# Saving to file
filename = 'X_train_file.csv' # Specify the file name

# Writing to csv file
with open(filename, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    
    # Writing the data
    csvwriter.writerows(X_train)
print(f"Data saved to {filename}")

#creating X_test
X_test=[] 
for item in paths_test:
    pdb_file = item[0] 
    sdf_file = item[1]
    #print(pdb_file)
    #print(sdf_file)  

    pro_names_cord  , protein_atoms = extract_atom_positions_pdb (pdb_file)# getting cords and store in protein_atoms
    lig_names_cord,ligand_atoms=extract_atom_positions_sdf(sdf_file) #getting cords and store in ligad_atoms
    cutoff_distance = 9.1  # Defining a suitable cutoff distance

    close_pro,lig = siever_fxn(pro_names_cord,lig_names_cord,cutoff_distance)
    features = features_extractor(close_pro,lig)
    X_test.append(features)
        
    print(len(features))
    print(len(X_test))
#print(X_test)

#saving to a file
file2name = 'X_test_file.csv'# Specify the file name

# Writing to csv file
with open(file2name, 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    
    # Writing the data
    csvwriter.writerows(X_test)

print(f"Data saved to {file2name}")

"""