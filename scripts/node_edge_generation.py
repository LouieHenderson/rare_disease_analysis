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

def genProtFuncEdge(funcod, protcod):
    #Generate Protein <> Functional edges
    outcsv = [[":START_ID",":END_ID", ":TYPE"]]
    count = 0
    for i in funcod:
        if funcod[i] in protcod:
            outcsv.append([i, protcod[funcod[i]], "is_feature"])

    return outcsv

def genDisProtEdge(disdic, protcod, discod, ensdic):
    #Generate Disease <> Functional edges
    outcsv = [[":START_ID",":END_ID", ":TYPE", "score"]]
    for i in disdic:
        for j in disdic[i]:
            if j["targetId"] in ensdic:
                outcsv.append([protcod[ensdic[j["targetId"]]], discod[i], "direct_evidence", float(j["score"])])
    return outcsv