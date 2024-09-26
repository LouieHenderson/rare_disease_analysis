# Rare disease analysis

Rare diseases are conditions that affect a  percentage of the population, typically fewer than 200,000 people in the U.S. or 1 in 2,000 in Europe. This imprecision in the definition makes determining the precise number of "rare diseases" , however estimates and curation efforts cite approximately 7000 known diseases.

The low incidence within the general population of these rare diseases make them economically unappealing for traditional drug discovery processes. This unmet need may be overcome through the repurposing of existing drugs to modulate and improve patient outcomes in understudied rare disease areas.


- https://health.ec.europa.eu/medicinal-products/orphan-medicinal-products_en    
- https://www.sciencedirect.com/science/article/pii/S0012369218300643


## Goal/s 

1. Generate graph database of rare disease ontologies linked to protein intel
2. Identify / integrate relevant datasets into graph db
3. Perform network analysis to identify drugs to repurpose for rare disease targets

## To-do

- Identify core datasets
- Parse datasets  
- Initialise graph
- Identify potential tractable targets / diseases 

## Datasets

### Protein intel

- Uniprot 
- OpenTargets
    - "targets"

### Disease intel / ontology

- Orphanet (ORPHA)
    - https://www.orphadata.com/orphanet-scientific-knowledge-files  
    - A rare disease nomenclature / ontology resource. Contains curated files crossreferencing a subset of rare diseases across different ontologies. 
- Disease Ontology (DOID)
- OpenTargets
    - "associationByOverallDirect"
    - "associationByOverallIndirect"