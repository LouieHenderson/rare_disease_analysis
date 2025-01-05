import subprocess
import  pandas as pd
import csv
import os
from neo4j import GraphDatabase

def nodeToCSVprep(genenode, diseasenode):
    #Generate out lists for .csv out
    geneout = [["id:ID", ":LABEL", "name", "ensemble_id"]]
    for i in genenode:
        geneout.append([genenode[i]["ID"], "Gene", i, genenode[i]["Ensemble"]])
        
    disout = [["id:ID", ":LABEL", "name", "uri"]]
    for i in diseasenode:
        disout.append([diseasenode[i]["ID"], "Disease", i, diseasenode[i]["URI"]])
        
    return geneout, disout

def writecsv(infile, outfile):
    #Write out csv
    with open(outfile, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=";")
        writer.writerows(infile)

def run_query(query, driver):
    #Submit query to graph database, retrieve output as df
    with driver.session() as session:
        result = session.run(query)
        rows = [record.data() for record in result]
        df = pd.DataFrame(rows)
        return df
    
def startGraph(passloc, directory):
    #Password
    with open(passloc, "r") as f:
        password = f.read()
    
    curdir = os.path.abspath(directory)
    #Commands for stopping / generating / starting neo4j
    cmd0 = ["sudo", "-S", "neo4j", "stop"]
    cmd1 = [ \
            "sudo", "-S", "neo4j-admin", "database", "import", "full", "neo4j", \
            "--nodes="+ curdir+"/protnodes.csv", \
            "--nodes="+ curdir+"/disnodes.csv", \
            "--nodes="+ curdir+"/funcnodes.csv", \
            "--relationships="+ curdir+"/profuncedge.csv", \
            "--relationships="+ curdir+"/disprotedge.csv", \
            "--delimiter=;", \
            "--array-delimiter=|", \
            "--overwrite-destination", \
            "--verbose"
            ]
    cmd2 = ["sudo", "-S", "neo4j", "start", "--verbose"]
    
    #Stop server
    result = subprocess.run(cmd0, input=password+"\n", capture_output=True, text=True)
    if result.returncode == 0:
        print("Graph stopped")
    else:
        print(f"Stop failed with error:\n{result.stderr}")
    
    #Generate graph
    result = subprocess.run(cmd1, input=password+"\n", capture_output=True, text=True)
    if result.returncode == 0:
        print("Database import successful!")
    else:
        print(f"Import failed with error:\n{result.stderr}")
    
    #Start graph
    result = subprocess.run(cmd2, input=password+"\n", capture_output=True, text=True)
    if result.returncode == 0:
        print("Graph started")
    else:
        print(f"Start failed with error:\n{result.stderr}")