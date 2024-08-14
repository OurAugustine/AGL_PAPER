def get_data_dict(path: str) -> dict:
    # Initialize dicts to store data
    data_dict={}
    # Open the file and read line by line, skipping the first 5 lines
    with open(path, 'r') as file:
        for i, line in enumerate(file):
            if i < 5:
                continue  # Skip the first 5 lines
            # Parse each line to extract information
            fields = line.strip().split()  # Split line by whitespace
            
            # Store information in a dictionary
            data_dict.update({fields[0]: fields[3]})
    return data_dict
def sep_core_non_dup(ref_path:str,core_path:str)->list:
    """"this function takes in file path of the refined proteins  and core the return alist of proteing in the
    refined file and onther list o proteins that are not repeated in the core as alist"""
    refined_proteins_dict=get_data_dict(ref_path)
    core_proteins_dict=get_data_dict(core_path)

    refined_proteins_names= set(list(refined_proteins_dict.keys()))
    core_proteins_names=set(list(core_proteins_dict.keys()))
    non_duplicates= refined_proteins_names-core_proteins_names
    #we will use the refined train and non duplicted testing

    return list(refined_proteins_names), list(non_duplicates),list(core_proteins_names)

def main():
    ref_path="datasets/pdbbind_v2007/v2007/INDEX.2007.refined.data"  
    core_path="datasets/pdbbind_v2007/v2007/INDEX.2007.core.data" 
    all_folders,non_duplicates ,core_folder= sep_core_non_dup(ref_path,core_path)
    print('the count of all folders :',len(all_folders))
    print('count of repeated fonders:',len(non_duplicates))
    print('the leght of core folder:',len(core_folder))

    print('now printing the none repeated proteins')
    #print(non_duplicates)
    print(len(non_duplicates))
    #extrating none repeating items only from refined dict of proteing
    refined_proteins=get_data_dict(ref_path)
    filtered_proteins = {key: value for key, value in refined_proteins.items() if key in non_duplicates }
    #print(filtered_proteins)
    print(len(filtered_proteins))#just to ensure we have  all proteins we need
    dict_all=get_data_dict(ref_path)
    dict_core=get_data_dict(core_path)
    print(len(dict_all))
    print(dict_core)

if __name__ == "__main__":
    main()

