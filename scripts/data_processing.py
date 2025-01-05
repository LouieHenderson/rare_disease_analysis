from pyteomics import fasta
import pandas as pd
import numpy as np
import re
import time
from datetime import datetime
import json
import os
import subprocess

def genUniprot(unidata, genein):
    #Generate Uniprot <> Gene dictionary
    with open(genein) as f:
        read = f.readlines()
    unigene = {}
    for i in read:
        splitr = i.split("\t")
        if splitr[1] == "Gene_Name" or splitr[1] == "Gene_Synonym":
            unigene.setdefault(splitr[0], [])
            unigene[splitr[0]].append(splitr[2].strip("\n"))

    #Generate df of human proteome & dic of gene names to uniprot ids
    unidic    = {}
    fastain   = fasta.UniProt(unidata)
    count = 0
    for uni in fastain:
        count += 1
        unid = uni[0]['id']
        name = uni[0]['name']
        try:
            gene = unigene[unid][0]
        except:
            gene = uni[0]['gene_id']
        seq  = uni[1]
        unidic.setdefault(unid, {"Uniprot": unid, "Gene": gene, "Name": name, "Seq": seq})

    unidf = pd.DataFrame.from_dict(unidic, orient='index')

    return unidf

def parseOTtarget(data_dir, df):
    #Parse OpenTargets "targets" intel, return df
    #Add new intel cols to df
    intelcol = {"Ensembl_gene": 'id', \
                "Ensembl_trans": 'transcriptIds', \
                "Description": "functionDescriptions", \
                "Subcell_loc": "subcellularLocations", \
                "DB_entries": "dbXrefs", \
                "Bio_path": "pathways", \
                "Tx_approaches": "tractability"}
    for col in intelcol:
        df[col] = [list() for x in range(len(df.index))]
    
    data = []
    #Read in multi part .json to data list 
    for file in os.listdir(data_dir):
        path = data_dir + "/" + str(file).replace("._", "")
        with open(path, "r") as f:
            for line in f:
                data.append(json.loads(line))
                
    #Retrieve intel from protein coding entries, annotate df
    ensemtouni = {}
    for entry in data:
        try:
            entid = entry['proteinIds'][0]['id']
        except:
            entid = ""

        if entid in df.index:
            for j in intelcol:
               #if intelcol[j] in entry and j != "Bio_path":
                if intelcol[j] in entry and j not in ["Bio_path", "Subcell_loc"]:
                    df.loc[entid, j].append(entry[intelcol[j]])
                #Simplify biological pathways for easier querying
                elif intelcol[j] in entry and j == "Bio_path":
                    for path in entry[intelcol[j]]:
                        df.loc[entid, j].append([path['pathwayId'], path['pathway']])
                elif j == "Ensembl_gene":
                    ensemtouni.setdefault(entry[intelcol[j]], entid) 
                elif intelcol[j] in entry and j == "Subcell_loc":
                    allsubloc = []
                    if len(entry[intelcol[j]]) > 0:
                        for sub in entry[intelcol[j]]:
                            allsubloc.append(sub["location"].lower())   
                    df.loc[entid, j].append(allsubloc[0:])
                   # df.loc[entid, j] = allsubloc
                    
    
    return df, ensemtouni

def parseOTdisease(dataloc):
    #Parse OT disease <> protein evidence
    data = {}
    for file in os.listdir(dataloc):
        if ".json" in str(file):
            with open(dataloc + "/" + str(file).replace("._", "")) as f:
                for line in f.readlines():
                    linedict = json.loads(line)
                    data.setdefault(linedict["diseaseId"], [linedict])
                    data[linedict["diseaseId"]].append(linedict)
    return data

def findDisease(database):
    #Parse high level OT disease intel
    data = {}
    count = 0
    for file in os.listdir(database):
        if ".json" in str(file):
            count += 1
            try:
                with open(database + "/" + str(file)) as f:
                    for line in f.readlines():
                        linedict = json.loads(line)
                        data.setdefault(linedict["id"], linedict)
            except:
                pass

    return data

def funcGen(unisprot):
    #Return a dictionary of Uniprot to function annotations
    outdic = {}
    with open(unisprot, "r") as f:
        infile = f.read().replace(";", " ")
        splitin= infile.split("//\n")
    #Iterate per entry, split block, retrieve Uniprot + individual functions
    for uni in splitin:
        count = 0
        allfunc = [line for line in uni.split("\n") if line.startswith('FT')]
        allac = [line for line in uni.split("\n") if line.startswith('AC')]
        try:
            uniprot = re.search(r"\bAC\s+([A-Za-z0-9]+)", uni).group().split()[1]
        except:
            uniprot = ""
        funcdic = {}
        func = {}
        curannot = ""
        #Iterate over functions, append to out dic
        for i in allfunc:
            match = re.match(r"^\S+(\s+)\S+", i)
            if len(match.group(1)) == 3:
                funcdic.setdefault(count, func)
                curannot = ""
                spliti = i.split()
                count += 1
                nums = re.findall(r'\d+', spliti[-1])
                nums = [int(num) for num in nums]
                func = {"type": spliti[1].replace(";", " "), "range": nums, "note": "", "evidence": "", "id": ""}
            else:
                if re.search(r"(?<=\/)\w+=+", i) != None:
                    search = re.search(r"(?<=\/)\w+=+", i)
                    curannot = search.group()[:-1]
                    func[curannot] = [i.split("=")[1].replace("\"", "")] 
                    spliti = i.split("/id=")
                else:
                    func[curannot] = func[curannot][0]+" ".join(i.split()[1:])
        outdic.setdefault(uniprot, funcdic)
    return outdic

def parseOrpha(datain):
    ##TO-DO
        #Alter output to df with cols = ontologies    
    #Read in Orphanet json, return dataframe of Orpha code / Name / other database crossreferences
    with open(datain) as f:
        dataj = json.load(f)
    outdic = {}
    for i in dataj["JDBOR"][0]["DisorderList"][0]["Disorder"]:
        orpha = i["OrphaCode"]
        nam   = i["Name"][0]["label"]
        xrefs = {}
        try:
            for j in i["ExternalReferenceList"][0]['ExternalReference']:
                xrefs.setdefault(j['Source'], {"ID": j['id'], "Mapping": j['DisorderMappingRelation'][0]["Name"][0]["label"]})
        except:
            pass
        outdic.setdefault(orpha, {"Name": nam, "Xref": xrefs})

    outdf = pd.DataFrame.from_dict(outdic)

    return outdf.transpose() 

def assocUniprotDisease(datadic, ensgenetouni, diseasexref):
    #Generate dictionary linking Uniprot ID to disease associations + scores
    uniprotdisease = {}
    for i in datadic:
        for ent in datadic[i]:
            try:
                Orpha   = int(ent['diseaseId'].split("_")[1])
                Score   = round(ent['score'], 3)
                Uniprot = ensgenetouni[ent['targetId']]
                Disnam  = diseasexref.loc[str(Orpha)]["Name"]
                uniprotdisease.setdefault(Uniprot, {Orpha: {"Disease": Disnam, "Score": Score}})
                uniprotdisease[Uniprot][Orpha] = {"Disease": Disnam, "Score": Score}
            except:
                pass
    return uniprotdisease