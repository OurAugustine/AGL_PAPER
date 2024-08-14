Below I will write the work of each script file 
In the first file called Proterinfunctions.py 
The ain of Proteinfunctions.py is to open file of refined protein and file of core proteins, then select those protein which are in the the refined bun not in the core and name 
this list a non duplicate , this is then used for testing, the whole group for refined protein is put in to a list names refine ,which will be used for training

For the second File called Atom_position.py 
The file has two function  one that get atom cordinates from a protein file sent to it , and one that get atom cordinates from a ligand file sent to it.

The read_file.py ,Whose main aim it to eventually return a list of paths to each and every ligand and protein , given the main foler and protein names
The file as two function  
    a) find files :this function opens a folder and then return file path to files ending with .sdf and protein.pdb as a tuppleready to be appended in to a list
    b) get_path_list : this fxn takes in,path to main folder, list of folders in the main folder,extenstin of files you need to search in it, 
    then return complete path to all the files ending with the string  you have provided

The File ADJ-LAP-MAT.py :Main aim is to create the Weighted adjacency , degree and laplacian matrix in that order , then calculate eigen values and eigen vectors
