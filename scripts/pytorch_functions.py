
class SequenceClassifier(nn.Module):
    def __init__(self, embedding_dim, hidden_dim, num_classes):
        super(SequenceClassifier, self).__init__()
        self.rnn = nn.LSTM(embedding_dim, hidden_dim, batch_first=True)
        self.fc = nn.Linear(hidden_dim, num_classes)

    def forward(self, x):
        #print(x)
        _, (hidden, _) = self.rnn(x)  # hidden: [1, batch_size, hidden_dim]
        hidden = hidden[-1]#.squeeze(0)    # Remove sequence dimension
        output = self.fc(hidden)
        return output

def EmbeddingGeneration(inputlist, outdir, model):
    embdic = {}
    baduniprot = []
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    with torch.no_grad():
        #Load model to gpu
        if torch.cuda.is_available():
            model = model.cuda()
            
        batch_labels = None 
        batch_strs   = None
        batch_tokens = None
        
        start = time.time()
        #Tokenise + pad sequences
        batch_labels, batch_strs, batch_tokens = batch_converter(inputlist)
        #Dictionary for matching Uniprot <> Token
        comdic = dict(zip(batch_labels, batch_tokens))
    
        truestart = time.time()
        #Iterate over chunks of seqs
        for ent in range(len(batch_tokens)):
            start = time.time()
            #Clear GPU memory
            out          = None
            batch_subset = None
            gc.collect()
            torch.cuda.empty_cache()
            if os.path.exists(outdir+batch_labels[ent]+".pt") == False:
                try:
                    batch_subset = batch_tokens[ent].to(device="cuda", non_blocking=True)
                    out = model(batch_subset.unsqueeze(0), repr_layers=[33], return_contacts=False)["representations"][33]
                    out = out.cpu()
                    torch.save(out, outdir+batch_labels[ent]+".pt")
                except:
                    print(f"\nIssue! = {chunk, batch_labels[ent]}")
                    print(f"Pre - Allocated memory: {torch.cuda.memory_allocated() / 1024**2:.2f} MB")
                    print(f"Pre - Cached memory: {torch.cuda.memory_reserved() / 1024**2:.2f} MB")
                    baduniprot.append(batch_labels[ent])
            else:
                print(f"Exists - {batch_labels[ent]+".pt"}")
            print(f"Time taken (Embedding (batch {ent})) = {round(time.time()-start, 2)} seconds \t| Total ({round(time.time()-truestart, 2)})")

    return embdic, baduniprot

def loadEmbed(embdir, labs, maindf, limit):
    embdic = {}
    limit = limit
    count = 0
    for path, dirs, files in os.walk(embdir):
        for file in files:
            if file[-2:] == "pt" and count < limit:
                count += 1
                gc.collect()
                torch.cuda.empty_cache()
                lab = labs[maindf.loc[maindf["Uniprot"] == file[:-3]]["Location"].values[0]]
                print(file, f"\t|\t{count}")
                embdic.setdefault(file[:-3],  [lab, torch.load(embdir+file)])
                
    return embdic

