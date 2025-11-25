def genProtNode(df):
    #Generate Uniprot Protein node
    count = 0
    outcsv = [["id:ID", ":LABEL", "uniprot", "gene", "name", "seq", "ens_gene", "ens_trans", "description", "db_entries", "bio_path", "tx_approach", "subcellular_location"]]
    nodecodec = {}
    for index, row in df.iterrows():
        count += 1
        outcsv.append(["P"+str(count), "Protein", row["Uniprot"], row['Gene'], row['Name'], row['Seq'], row['Ensembl_gene'], row['Ensembl_trans'],\
                        row['Description'], row['DB_entries'], row['Bio_path'], row['Tx_approaches'], row["Subcell_loc"]])
        nodecodec.setdefault(index, "P"+str(count))
    return outcsv, nodecodec

def genFuncNode(funcdic):
    #Generate Uniprot Function node
    count = 0
    outcsv = [["id:ID", ":LABEL", "type",  "range", "alt_id", "evidence", "note"]]
    nodecodec = {}
    for uni in funcdic:
        for func in funcdic[uni]: 
            if len(funcdic[uni][func]) > 0:
                count += 1
                outcsv.append(["F"+str(count), "Function", funcdic[uni][func]["type"], funcdic[uni][func]["range"], funcdic[uni][func]["id"], funcdic[uni][func]["evidence"], funcdic[uni][func]["note"]])
                nodecodec.setdefault("F"+str(count), uni)
    return outcsv, nodecodec

def genDisNode(disdic, disedge):
    #Generate OpenTarget disease node
    count = 0
    outcsv = [["id:ID", ":LABEL", "id_ont", "name", "description"]]
    nodecodec = {}
    for dis in disdic:
        if dis in disedge:
            try:
                desc = disdic[dis]["description"]
            except:
                desc = ""
            count += 1
            outcsv.append(["D"+str(count), "Disease", dis, disdic[dis]["name"], desc.replace("\n", " ")])
            nodecodec.setdefault(dis, "D"+str(count))
    return outcsv, nodecodec

def genInDisNode(indisdic, indisedge):
    #Generate OpenTarget disease node
    count = 0
    outcsv = [["id:ID", ":LABEL", "id_ont", "name", "description"]]
    nodecodec = {}
    for dis in indisdic:
        if dis in indisedge:
            try:
                desc = indisdic[dis]["description"]
            except:
                desc = ""
            count += 1
            outcsv.append(["ID"+str(count), "Indirect_Disease", dis, indisdic[dis]["name"], desc.replace("\n", " ")])
            nodecodec.setdefault(dis, "ID"+str(count))
    return outcsv, nodecodec

def genProtFuncEdge(funcod, protcod):
    #Generate Protein <> Functional edges
    outcsv = [[":START_ID",":END_ID", ":TYPE"]]
    count = 0
    for i in funcod:
        if funcod[i] in protcod:
            outcsv.append([i, protcod[funcod[i]], "is_feature"])

    return outcsv

def genDisProtEdge(disdic, protcod, discod, ensdic, mainlable):
    #Generate Disease <> Functional edges
    outcsv = [[":START_ID",":END_ID", ":TYPE", "score"]]
    for i in disdic:
        for j in disdic[i]:
            if j["targetId"] in ensdic:
                outcsv.append([protcod[ensdic[j["targetId"]]], discod[i], mainlable, float(j["score"])])
    return outcsv

def genProtSubEdge(protcod, subcelcod, df):
    #Generate Protein <> Subcellular location edges
    outcsv = [[":START_ID",":END_ID", ":TYPE"]]
    for sub in subcelcod:
        subdf = df[df["Subcell_loc"].apply(lambda x: sub in x[0] if x else False)]
        for index, row in subdf.iterrows():
            outcsv.append([protcod[row["Uniprot"]], subcelcod[sub], "is_location"])
    return outcsv

def genSubcellNode(df):
    #Subcellular localisation annotations for broad categories
    uniloc = {"Intracellular": [], "Extracellular": [], "Membrane": []}

    #Identify unique locations
    indiv_loc = {}
    for loc_list in df["Subcell_loc"].values:
        for ent in loc_list:
            for loc in ent:
                if loc not in indiv_loc:
                    indiv_loc.setdefault(loc, "")
                    
    #Add broad subcellular locations based on keywords
    for loc in indiv_loc:
        if any( i in loc.lower().split() for i in ["secreted", "secret", "extracellular"]):
            indiv_loc[loc] = "extracellular"
        else:
            if any( i in loc.lower().split() for i in ["membrane", "wall", "lipid"]):
                indiv_loc[loc] = "membrane"
            else:
                indiv_loc[loc] = "intracellular"
    
    #Generate Subcellular Location node
    count = 0
    outcsv = [["id:ID", ":LABEL", "location", "intra_mem_extra"]]
    nodecodec = {}
    for loc in indiv_loc:
        count += 1
        outcsv.append(["S"+str(count), "Subcellular", loc, indiv_loc[loc]])
        nodecodec.setdefault(loc, "S"+str(count))

    return outcsv, nodecodec