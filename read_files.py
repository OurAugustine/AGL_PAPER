from Proteinfunctions import sep_core_non_dup
import os

def find_files(directory:str, pdb_extension:str,sdf_extension:str)-> tuple:
    """this function opens a folder and then return file path to files ending with .sdf and protein.pdb as a tupple
    ready to be appended in to a list """
    
    pdb_path=''
    sdf_path=''
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(sdf_extension):
                sdf_path= directory + "/"+file
            elif file.endswith(pdb_extension):
                pdb_path= directory +"/"+ file
    
    return (pdb_path,sdf_path) #file_paths




def get_path_list(main_folder:str,folder_list:list,end1:str,end2:str)->list:
    '''this fxn takes in,path to main folder, list of folders in the main folder,extenstin of files you need to search 
    in it, then return complete path to all the files ending with the string  you have provided'''
    
    all_path_list=[]
    for a_folder in folder_list:
        path = main_folder + "/" + a_folder #create a path'
        filepaths=find_files(path,end1,end2)# use the path to comple path to our needed files
        all_path_list.append(filepaths) #add the two paths to our list
    return all_path_list
    
# Usage
def main():
    ref_path="datasets/pdbbind_v2007/v2007/INDEX.2007.refined.data"  
    core_path="datasets/pdbbind_v2007/v2007/INDEX.2007.core.data"
    all_folder_list,non_dup_folder_list=sep_core_non_dup(ref_path,core_path)
    #returns alist of :all refined atoms,refined but non duplicated atoms
    main_folder = "datasets/pdbbind_v2007/v2007" 
    extention1 = 'n.pdb'  # replace with your file extension
    extention2= '.sdf'
    paths= get_path_list(main_folder,non_dup_folder_list,'n.pdb','.sdf')
    print(len(paths))
    print(paths[0:5])

if __name__ == "__main__":
    main()