from multiprocessing import Process

def wrapper(id,subset_size,protein_length,num_procs):
    # generate file names for work
    temp_file = "temp_file_" + str(id)
    out_file = 'outfile_' + str(id)
    #find the number of folders to look at
    start_step = id
    #quick check to make sure we grab every folder, last process might have more folders than the rest BUT later folders are smaller so that's not bad
    #We could try and add some math to distribute the load more evenly, like give each process every ith valued folder
        
    folders_to_analyze = list(range(start_step,protein_length,num_procs))
    if id == num_procs-1:
        folders_to_analyze.append(protein_length)
    
    #this is our final list with what we should need
    values = [temp_file,out_file,folders_to_analyze]
    print(values)
    #main(value)
    

if __name__ == '__main__':
    # Just boot up however many processes we want and have them start the function
    # the math in the wrapper covers generating names for all the temp and output files as well as how many folders each process needs to examine
    # all the values we should have to touch are the number of processes and the protein protein_length
    # as well as any paths used anywhere
    num_procs = 20
    protein_length = 101
    for i in range (0,num_procs):
        p = Process(target=wrapper, args=(i,int(protein_length/num_procs),protein_length,num_procs))
        p.start()
        p.join()
