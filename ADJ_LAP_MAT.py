import numpy as np
from atom_position import extract_atom_positions_pdb # return name, cords
from atom_position import extract_atom_positions_sdf #returns name, cords
from atom_position import siever_fxn

# 1. AGL paper provides the algebraic graph laplacian matrix
# 2. Calculate all the features for the PDBBind 2007 dataset
# 3. Use Gradient Boosting Regressor to predict the binding affinity
#    - Use the refined dataset as the training set
#    - Use the code dataset as the test set
#    - Make sure that the intersection of the two datasets is empty
#      If not, remove the samples in the refined dataset from the code dataset

#coombined adj and lap matrxies
print("adjacency , degree and laplacian matrix from combbined data set of sdf and pdb")
def create_combined_adjacency_matrix(protein_atoms: list[tuple[float, float, float]],
                                     ligand_atoms :list[tuple[float, float, float]], cutoff: float) -> np.ndarray:
    
    total_atoms = len(protein_atoms) + len(ligand_atoms)
    adjacency_matrix = np.zeros((total_atoms, total_atoms), dtype=float)
       
    for protein_atom in protein_atoms:
        for ligand_atom in ligand_atoms:
            distance = np.linalg.norm(np.array(protein_atom) - np.array(ligand_atom))
            if distance <= cutoff:
                protein_index = protein_atoms.index(protein_atom)
                ligand_index = len(protein_atoms) + ligand_atoms.index(ligand_atom)
                adjacency_matrix[protein_index, ligand_index] = np.exp(-distance)
                adjacency_matrix[ligand_index, protein_index] = np.exp(-distance)
    
    return adjacency_matrix

def create_laplacian_matrix(adjacency_matrix: np.ndarray) -> np.ndarray:
    degree_matrix=np.diag(np.sum(adjacency_matrix,axis=1))
    laplacian_matrix = degree_matrix - adjacency_matrix
    return laplacian_matrix


def fetch_eig_stat(matrix:np.ndarray)->np.ndarray:
    """this function takes in a matrix statistics of its eigien value as 1)1st non zero ,2)sum ,3)count of zeros ,4)count of non zeros, 
    5)mean of positve eig val,6)standard division of +ve eig vals,7)min ,8) max , and 9) dot product of eigen vals"""
    eigenvalues= np.linalg.eigvals(matrix)
    pos_eig_vals=[val for val in eigenvalues if val >=0]
    none_zeros = [val for val in eigenvalues if val != 0]
    #min_non_zero= min(none_zeros)
    #max_non_zero=max(none_zeros)
    first_non_zero=next((x for x in eigenvalues if x != 0), 0)
    zeros_count = list(eigenvalues).count(0)
    non_zero_count = len(eigenvalues) - zeros_count
    sumed=sum(pos_eig_vals)

    minmum_val=0
    max_val = 0
    if len(eigenvalues) == 0:
        minmum_val=0
        max_val=0
    else:
        minmum_val=min(eigenvalues)
        max_val = max(eigenvalues)

    meaned=0
    std=0
    if len(pos_eig_vals) == 0:
        meaned=0
        std=0
        #minmum_val=0
    else:
        meaned=np.mean(pos_eig_vals)
        std= np.std(pos_eig_vals)
        #minmum_val=min(eigenvalues)

    dot_prod=np.dot(eigenvalues,eigenvalues)
    #list1=[first_non_zero,sumed,zeros_count,non_zero_count,meaned,std,minmum_val,max_val,dot_prod]   
    return first_non_zero,sumed,zeros_count,non_zero_count,meaned,std,minmum_val,max_val,dot_prod #list1


def features_extractor(pro_atoms:list,lig_atoms:list)->list:
    """"this fnct takes in 2 lists both in the format list=[(string(float,float, float)),....]], it the creates 36 small matrixes for all interaction
      between protein elements and ligand elements of intres , there after  calculate eigen values for each matrix and retrun there statistic """    
    features = []
    pro_elements=['O','S','N','C']
    for i_type in pro_elements:
        lig_elements = ["H","C","N","O","S","P","F","Cl","Br"]
        for j_type in lig_elements:
            filtered_pro_atoms = [b for a,b in pro_atoms if a == i_type]
            filtered_lig_atoms = [b for a,b in lig_atoms if a == j_type]
            total_atoms = len(filtered_pro_atoms) + len(filtered_lig_atoms)
            adjacency_matrix = np.zeros((total_atoms, total_atoms), dtype=float)

        
            for protein_atom in filtered_pro_atoms:
                for ligand_atom in filtered_lig_atoms:
                    distance = np.linalg.norm(np.array(protein_atom) - np.array(ligand_atom))
                    protein_index = filtered_pro_atoms.index(protein_atom)
                    ligand_index = len(filtered_pro_atoms) + filtered_lig_atoms.index(ligand_atom)
                    adjacency_matrix[protein_index, ligand_index] = np.exp(-distance)
                    adjacency_matrix[ligand_index, protein_index] = np.exp(-distance)
        
            laplacian_matrix= create_laplacian_matrix(adjacency_matrix)
            #eigenvalues,eigenvectors=np.linalg.eig(laplacian_matrix)
            #features.append(sum(eigenvalues))
            list1=fetch_eig_stat(laplacian_matrix)
            for item in list1:
                    features.append(item)            
    return features

def main():
    # Usage example
    pdb_file = "datasets/pdbbind_v2007/v2007/1a0q/1a0q_protein.pdb" 
    sdf_file = "datasets/pdbbind_v2007/v2007/1a0q/1a0q_ligand.sdf"   

    pro_atoms, protein_atoms = extract_atom_positions_pdb (pdb_file)# getting cords and store in protein_atoms
    #ligand_atoms, _, _ = parse_sdf(sdf_file)#getting cords and store in ligad_atoms
    lig_atoms,ligand_atoms=extract_atom_positions_sdf(sdf_file)#getting cords and store in ligad_atoms
    cutoff_distance = 4.0  # Defining a suitable cutoff distance

    adjacency_matrix = create_combined_adjacency_matrix(protein_atoms, ligand_atoms, cutoff_distance)
    laplacian_matrix = create_laplacian_matrix(adjacency_matrix)

    print("Adjacency Matrix:")
    print(adjacency_matrix)
    print("\nLaplacian Matrix:")
    print(laplacian_matrix)
    print('\n the eigen values statiscits')
    stats=fetch_eig_stat(laplacian_matrix)
    print(stats)
    
    print('working withs smaller matrix')
    close_pro,lig=siever_fxn(pro_atoms,lig_atoms,cut_off=9.5)
    features =features_extractor(close_pro,lig)
    print(features)
    print(len(features))
if __name__ == "__main__":
    main()