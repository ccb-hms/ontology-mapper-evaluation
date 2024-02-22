import os
import time
import bioregistry
import pandas as pd
import numpy as np
import text2term
from text2term import Mapper
from tqdm import tqdm


__version__ = "0.1.2"

# URL to EFO ontology version used
EFO_URL = "http://www.ebi.ac.uk/efo/releases/v3.62.0/efo.owl"

# Output folder for the mappings computed using text2term and those extracted from GWAS Catalog
OUTPUT_FOLDER = "output"

# Headers of the GWAS Catalog metadata table
STUDY_ID_COLUMN = "STUDY.ACCESSION"
TRAIT_COLUMN = "DISEASE.TRAIT"
MAPPED_TRAIT_COLUMN = "MAPPED_TRAIT"
MAPPED_TRAIT_IRI_COLUMN = "MAPPED_TRAIT_URI"
MAPPED_TRAIT_CURIE_COLUMN = "MAPPED_TRAIT_CURIE"

# Headers of the text2term mappings table
T2T_INPUT_TERM_COL = "Source Term"
T2T_INPUT_TERM_ID_COL = "Source Term ID"
T2T_MAPPED_TERM_CURIE_COL = "Mapped Term CURIE"


def compute_text2term_mappings(metadata_df, source_term_col=TRAIT_COLUMN, source_term_id_col=STUDY_ID_COLUMN):
    source_terms = metadata_df[source_term_col].tolist()
    source_term_ids = metadata_df[source_term_id_col].tolist()
    return text2term.map_terms(source_terms=source_terms, source_terms_ids=source_term_ids,
                               target_ontology=EFO_URL, excl_deprecated=True, save_graphs=False,
                               max_mappings=1, min_score=0.0, save_mappings=True, mapper=Mapper.TFIDF,
                               output_file=os.path.join(OUTPUT_FOLDER, "mappings_t2t.csv"),
                               base_iris=("http://www.ebi.ac.uk/efo/", "http://purl.obolibrary.org/obo/MONDO",
                                          "http://purl.obolibrary.org/obo/HP", "http://www.orpha.net/ORDO",
                                          "http://purl.obolibrary.org/obo/DOID"))


def extract_gwascatalog_mappings(metadata_df):
    print("Extracting ontology mappings from the GWAS Catalog metadata...")
    mappings_list = []
    for _, row in metadata_df.iterrows():
        mapped_trait_uri = row[MAPPED_TRAIT_IRI_COLUMN]
        if mapped_trait_uri != "" and not pd.isna(mapped_trait_uri):
            if "," in mapped_trait_uri:
                iris_list = mapped_trait_uri.split(',')
            else:
                iris_list = [mapped_trait_uri]
            # Add a new row to the "mappings_df" for each IRI in the list
            for iri in iris_list:
                iri = iri.strip()
                mappings = {STUDY_ID_COLUMN: row[STUDY_ID_COLUMN],
                            TRAIT_COLUMN: row[TRAIT_COLUMN],
                            MAPPED_TRAIT_COLUMN: row[MAPPED_TRAIT_COLUMN],
                            MAPPED_TRAIT_IRI_COLUMN: iri,
                            MAPPED_TRAIT_CURIE_COLUMN: _get_curie_for_term(iri)}
                mappings_list.append(mappings)
    output_mappings_file = os.path.join(OUTPUT_FOLDER, "mappings_gwascatalog.tsv")
    mappings_df = pd.DataFrame(mappings_list)
    mappings_df.to_csv(output_mappings_file, sep="\t", index=False)
    print(f"...done (saved to {output_mappings_file})")
    return mappings_df


def _get_curie_for_term(term):
    if (not pd.isna(term)) and ("<" in term or "http" in term):
        term = term.replace("<", "")
        term = term.replace(">", "")
        if "," in term:
            tokens = term.split(",")
            curies = [_get_curie(token.strip()) for token in tokens]
            curies = ",".join(curies)
            return curies
        else:
            return _get_curie(term)
    return term


def _get_curie(term):
    curie = bioregistry.curie_from_iri(term)
    if curie is None:
        if "http://dbpedia.org" in term:
            return "DBR:" + term.rsplit('/', 1)[1]
        else:
            return term
    curie = curie.upper()
    if "OBO:" in curie:
        curie = curie.replace("OBO:", "obo:")
    if "NCBITAXON:" in curie:
        curie = curie.replace("NCBITAXON:", "NCBITaxon:")
    if "ORPHANET.ORDO" in curie:
        curie = curie.replace("ORPHANET.ORDO", "ORDO")
    return curie


def compare_mappings(t2t_mappings, gwascat_mappings, efo_df):
    data = []
    unique_accessions = t2t_mappings[T2T_INPUT_TERM_ID_COL].unique()

    for accession in tqdm(unique_accessions):
        t2t_subset = t2t_mappings[t2t_mappings[T2T_INPUT_TERM_ID_COL] == accession]
        gwascat_subset = gwascat_mappings[gwascat_mappings[STUDY_ID_COLUMN] == accession]

        input_trait = t2t_subset[T2T_INPUT_TERM_COL].iloc[0]
        t2t_traits = t2t_subset[T2T_MAPPED_TERM_CURIE_COL].unique()
        gwascat_traits = gwascat_subset[MAPPED_TRAIT_CURIE_COLUMN].unique()

        t2t_trait = t2t_traits[0]  # there is a single mapping for each input trait
        t2t_trait_label = t2t_subset["Mapped Term Label"].unique()
        if not isinstance(t2t_trait_label, str):
            t2t_trait_label = t2t_trait_label[0]

        # look up parents and children in the EFO class hierarchy
        trait_parents = efo_df.loc[efo_df['Subject'] == t2t_trait, 'Object'].unique()
        trait_children = efo_df.loc[efo_df['Object'] == t2t_trait, 'Subject'].unique()

        if np.any(gwascat_traits[:] == t2t_trait):
            data.append([accession, input_trait, t2t_trait, t2t_trait_label, t2t_trait, t2t_trait_label, 'Same'])
        elif any(t in trait_parents for t in gwascat_traits):
            gwascat_trait_label = gwascat_subset[MAPPED_TRAIT_COLUMN].iloc[0]
            data.append([accession, input_trait, t2t_trait, t2t_trait_label, gwascat_traits[0], gwascat_trait_label[0], 'More Specific'])
        elif any(t in trait_children for t in gwascat_traits):
            gwascat_trait_label = gwascat_subset[MAPPED_TRAIT_COLUMN].iloc[0]
            data.append([accession, input_trait, t2t_trait, t2t_trait_label, gwascat_traits[0], gwascat_trait_label[0], 'More Generic'])
        else:
            gwascat_trait = ""
            gwascat_trait_label = ""
            if len(gwascat_traits) > 0:
                gwascat_trait = gwascat_traits[0]
                gwascat_trait_label = gwascat_subset[MAPPED_TRAIT_COLUMN].iloc[0]
            data.append([accession, input_trait, t2t_trait, t2t_trait_label, gwascat_trait, gwascat_trait_label, 'Unrelated'])

    return pd.DataFrame(data, columns=[STUDY_ID_COLUMN, TRAIT_COLUMN, "TEXT2TERM.MAPPING", "TEXT2TERM.MAPPING.LABEL",
                                       "GWASCAT.MAPPING", "GWASCAT.MAPPING.LABEL", "CATEGORY"])


if __name__ == '__main__':
    # Create output directory if it does not exist
    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)

    # Load the GWAS Catalog metadata table
    gwascatalog_metadata = pd.read_csv(os.path.join("data", "gwascatalog_metadata.tsv"), sep="\t")

    # Extract ontology mappings from the GWAS Catalog metadata
    gwascatalog_mappings = extract_gwascatalog_mappings(metadata_df=gwascatalog_metadata)

    # Load the EFO ontology table contained all entailed SubClassOf relationships between terms in EFO
    efo_edges_df = pd.read_csv(os.path.join("data", "efo_entailed_edges.tsv"), sep="\t")

    # Compute text2term mappings from scratch or load mappings from file if they exist in the OUTPUT_FOLDER
    t2t_mappings_file = os.path.join(OUTPUT_FOLDER, "mappings_t2t.csv")
    if os.path.exists(t2t_mappings_file):
        text2term_mappings = pd.read_csv(t2t_mappings_file, skiprows=11, low_memory=False)
        print("Loading text2term mappings from file...")
    else:
        text2term_mappings = compute_text2term_mappings(metadata_df=gwascatalog_metadata)

    start = time.time()
    print("Comparing mappings...")
    results_df = compare_mappings(text2term_mappings, gwascatalog_mappings, efo_edges_df)
    output_file = os.path.join(OUTPUT_FOLDER, "mappings_comparison.tsv")
    results_df.to_csv(output_file, sep="\t", index=False)
    print(f"...done (comparison time: {format(time.time() - start, '.2f')} seconds. Saved results to: {output_file})")
