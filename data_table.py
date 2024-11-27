import pandas as pd

file = "1crn_total.txt"
longest_line = -1
pandas_data = []
headings = []

with open(file, "r") as f:
    lines = f.readlines()
    for line in lines:
        split_line = line.split(" ")
        if len(split_line) > longest_line:
            longest_line = len(split_line)
        split_line[-1] = split_line[-1].strip("\n")
        pandas_data.append(split_line)
        
        
    print(len(pandas_data))
    headings = ["Insertion Location 1", "Residue 1","Insertion Location 2", "Residue 2"]
    for i in range(0,longest_line-4):
        headings.append("Hbond distance " + str(i))
    print(len(headings))

df = pd.DataFrame(pandas_data,columns=headings)
df.to_csv("1crn_hbonds.csv", index=False)

           


