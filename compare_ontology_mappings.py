import os
import time
import bioregistry
import pandas as pd
import text2term
from text2term import Mapper
from tqdm import tqdm

__version__ = "0.3.0"

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

# Dictionary specifying for each term their high-level types/parents
MAPPING_TYPES = dict()


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


def extract_gwascatalog_mappings(metadata_df, entailed_edges_df):
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
                term_curie = _get_curie_for_term(iri)
                term_parents = entailed_edges_df.loc[entailed_edges_df['Subject'] == term_curie, 'Object'].unique()

                is_disease = _has_parent("EFO:0000408", term_parents)
                is_measurement = _has_parent("EFO:0001444", term_parents)
                is_phenotype = _has_parent("EFO:0000651", term_parents)
                is_amount = _has_parent("PATO:0000070", term_parents)
                is_process = _has_parent("BFO:0000015", term_parents)

                MAPPING_TYPES[term_curie] = {"IS_DISEASE": is_disease,
                                             "IS_MEASUREMENT": is_measurement,
                                             "IS_PHENOTYPE": is_phenotype,
                                             "IS_AMOUNT": is_amount,
                                             "IS_PROCESS": is_process}

                mappings = {STUDY_ID_COLUMN: row[STUDY_ID_COLUMN],
                            TRAIT_COLUMN: row[TRAIT_COLUMN],
                            MAPPED_TRAIT_COLUMN: row[MAPPED_TRAIT_COLUMN],
                            MAPPED_TRAIT_IRI_COLUMN: iri,
                            MAPPED_TRAIT_CURIE_COLUMN: term_curie,
                            "IS_DISEASE": is_disease,
                            "IS_MEASUREMENT": is_measurement,
                            "IS_PHENOTYPE": is_phenotype,
                            "IS_AMOUNT": is_amount,
                            "IS_PROCESS": is_process,
                            "PARENTS": term_parents}
                mappings_list.append(mappings)
    output_mappings_file = os.path.join(OUTPUT_FOLDER, "mappings_gwascatalog.tsv")
    mappings_df = pd.DataFrame(mappings_list)
    mappings_df.to_csv(output_mappings_file, sep="\t", index=False)
    print(f"...done (saved to {output_mappings_file})")
    return mappings_df


def _has_parent(parent, parents):
    if parent in parents:
        return "TRUE"
    return "FALSE"


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


def compare_mappings(t2t_mappings, gwascat_mappings, edges_df, entailed_edges_df):
    data = []
    unique_accessions = t2t_mappings[T2T_INPUT_TERM_ID_COL].unique()

    for accession in tqdm(unique_accessions):
        t2t_subset = t2t_mappings[t2t_mappings[T2T_INPUT_TERM_ID_COL] == accession]
        gwascat_subset = gwascat_mappings[gwascat_mappings[STUDY_ID_COLUMN] == accession]

        input_trait = t2t_subset[T2T_INPUT_TERM_COL].iloc[0]
        t2t_trait = t2t_subset[T2T_MAPPED_TERM_CURIE_COL].iloc[0]  # there is a single mapping for each input trait
        t2t_trait_label = t2t_subset["Mapped Term Label"].iloc[0]

        gwascat_traits = gwascat_subset[MAPPED_TRAIT_CURIE_COLUMN].unique()

        # get asserted EFO parents of both the text2term mapped term and all gwascatalog mapped terms
        t2t_trait_asserted_parents = edges_df.loc[edges_df['Subject'] == t2t_trait, 'Object'].unique()
        gwascat_traits_asserted_parents = []
        for gwascat_trait in gwascat_traits:
            gwascat_trait_parents = edges_df.loc[edges_df['Subject'] == gwascat_trait, 'Object'].unique()
            gwascat_traits_asserted_parents.extend(gwascat_trait_parents)

        # get entailed parents and children (of the text2term mapped term) in the EFO class hierarchy
        t2t_trait_parents = entailed_edges_df.loc[entailed_edges_df['Subject'] == t2t_trait, 'Object'].unique()
        t2t_trait_children = entailed_edges_df.loc[entailed_edges_df['Object'] == t2t_trait, 'Subject'].unique()

        if len(gwascat_traits) > 0:
            gwascat_trait = gwascat_traits[0]
            gwascat_trait_label = gwascat_subset[MAPPED_TRAIT_COLUMN].iloc[0]

            trait_term = MAPPING_TYPES[gwascat_trait]
            is_disease = trait_term["IS_DISEASE"]
            is_measurement = trait_term["IS_MEASUREMENT"]
            is_phenotype = trait_term["IS_PHENOTYPE"]
            is_amount = trait_term["IS_AMOUNT"]
            is_process = trait_term["IS_PROCESS"]

            if not isinstance(gwascat_trait_label, str):
                gwascat_trait_label = gwascat_trait_label[0]
            if t2t_trait in gwascat_traits:
                data.append([accession, input_trait, t2t_trait, t2t_trait_label, gwascat_trait, gwascat_trait_label,
                             'Same', is_disease, is_measurement, is_phenotype, is_amount, is_process])
            elif any(t in t2t_trait_parents for t in gwascat_traits):
                data.append([accession, input_trait, t2t_trait, t2t_trait_label, gwascat_trait, gwascat_trait_label,
                             'More Specific', is_disease, is_measurement, is_phenotype, is_amount, is_process])
            elif any(t in t2t_trait_children for t in gwascat_traits):
                data.append([accession, input_trait, t2t_trait, t2t_trait_label, gwascat_trait, gwascat_trait_label,
                             'More General', is_disease, is_measurement, is_phenotype, is_amount, is_process])
            elif any(t in t2t_trait_asserted_parents for t in gwascat_traits_asserted_parents):
                data.append([accession, input_trait, t2t_trait, t2t_trait_label, gwascat_trait, gwascat_trait_label,
                             'Sibling', is_disease, is_measurement, is_phenotype, is_amount, is_process])
            else:
                data.append([accession, input_trait, t2t_trait, t2t_trait_label, gwascat_trait, gwascat_trait_label,
                             'Unrelated', is_disease, is_measurement, is_phenotype, is_amount, is_process])

    return pd.DataFrame(data, columns=[STUDY_ID_COLUMN, TRAIT_COLUMN, "TEXT2TERM.MAPPING", "TEXT2TERM.MAPPING.LABEL",
                                       "GWASCAT.MAPPING", "GWASCAT.MAPPING.LABEL", "CATEGORY", "IS_DISEASE",
                                       "IS_MEASUREMENT", "IS_PHENOTYPE", "IS_AMOUNT", "IS_PROCESS"])


if __name__ == '__main__':
    # Create output directory if it does not exist
    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)

    # Load the GWAS Catalog metadata table
    metadata_file = os.path.join("data", "gwascatalog_metadata.tsv")
    print(f"Loading GWAS Catalog metadata from: {metadata_file}")
    gwascatalog_metadata = pd.read_csv(metadata_file, sep="\t")

    # Filter out the studies/rows that have been mapped to multiple EFO ontology terms
    gwascatalog_metadata = gwascatalog_metadata[
        ~gwascatalog_metadata[MAPPED_TRAIT_CURIE_COLUMN].astype(str).str.contains(',')]
    print(f"...metadata contains {gwascatalog_metadata.shape[0]} traits")

    # Load the EFO ontology table containing all asserted SubClassOf relationships between terms in EFO
    efo_edges_df = pd.read_csv(os.path.join("data", "efo_edges.tsv"), sep="\t")

    # Load the EFO ontology table containing all entailed SubClassOf relationships between terms in EFO
    efo_entailed_edges_df = pd.read_csv(os.path.join("data", "efo_entailed_edges.tsv"), sep="\t")

    # Extract ontology mappings from the GWAS Catalog metadata
    gwascatalog_mappings = extract_gwascatalog_mappings(metadata_df=gwascatalog_metadata,
                                                        entailed_edges_df=efo_entailed_edges_df)

    # Compute text2term mappings from scratch or load mappings from file if they exist in the OUTPUT_FOLDER
    t2t_mappings_file = os.path.join(OUTPUT_FOLDER, "mappings_t2t.csv")
    if os.path.exists(t2t_mappings_file):
        print(f"Loading text2term mappings from file ({t2t_mappings_file})...")
        text2term_mappings = pd.read_csv(t2t_mappings_file, skiprows=11, low_memory=False)
    else:
        text2term_mappings = compute_text2term_mappings(metadata_df=gwascatalog_metadata)

    start = time.time()
    print("Comparing mappings...")
    results_df = compare_mappings(t2t_mappings=text2term_mappings, gwascat_mappings=gwascatalog_mappings,
                                  edges_df=efo_edges_df, entailed_edges_df=efo_entailed_edges_df)
    output_file = os.path.join(OUTPUT_FOLDER, "mappings_comparison.tsv")
    results_df.to_csv(output_file, sep="\t", index=False)
    print(f"...done (comparison time: {format(time.time() - start, '.2f')} seconds. Saved results to: {output_file})")
