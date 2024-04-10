import os
import time
import logging
import bioregistry
import pandas as pd
import text2term
from text2term import Mapper, onto_utils
from tqdm import tqdm

__version__ = "0.7.0"

# URL to EFO ontology version used
EFO_URL = "http://www.ebi.ac.uk/efo/releases/v3.62.0/efo.owl"

# Output folder for the mappings computed using text2term and those extracted from the chosen benchmark dataset
OUTPUT_FOLDER = "output"

# Headers of the text2term mappings table
T2T_INPUT_TERM_COL = "Source Term"
T2T_INPUT_TERM_ID_COL = "Source Term ID"
T2T_MAPPED_TERM_COL = "Mapped Term Label"
T2T_MAPPED_TERM_IRI_COL = "Mapped Term IRI"
T2T_MAPPED_TERM_CURIE_COL = "Mapped Term CURIE"

# Dictionary specifying for each term their high-level types/parents
MAPPING_TYPES = dict()

LOG = onto_utils.get_logger("mapping.comparator", logging.INFO)


def compare(benchmark_dataset, benchmark_dataset_name, benchmark_mappings, source_term_col, source_term_id_col):
    # Compute text2term mappings from scratch or load mappings from file if they exist in the output folder
    t2t_mappings_file = os.path.join(OUTPUT_FOLDER, f"{benchmark_dataset_name}_t2t_mappings.csv")
    if os.path.exists(t2t_mappings_file):
        LOG.info(f"Loading text2term mappings from file ({t2t_mappings_file})...")
        text2term_mappings = pd.read_csv(t2t_mappings_file, skiprows=11, low_memory=False)
    else:
        text2term_mappings = compute_text2term_mappings(metadata_df=benchmark_dataset,
                                                        dataset_name=benchmark_dataset_name,
                                                        source_term_col=source_term_col,
                                                        source_term_id_col=source_term_id_col)
    start = time.time()
    LOG.info("Comparing mappings...")
    results_df = compare_mappings(t2t_mappings=text2term_mappings, benchmark_mappings=benchmark_mappings,
                                  edges_df=efo_edges_df, entailed_edges_df=efo_entailed_edges_df)
    output_file = os.path.join(OUTPUT_FOLDER, f"{benchmark_dataset_name}_results.tsv")
    results_df.to_csv(output_file, sep="\t", index=False)
    LOG.info(
        f"...done (comparison time: {format(time.time() - start, '.2f')} seconds. Saved results to: {output_file})")


def compute_text2term_mappings(metadata_df, dataset_name, source_term_col, source_term_id_col):
    source_terms = metadata_df[source_term_col].tolist()
    source_term_ids = metadata_df[source_term_id_col].tolist()
    return text2term.map_terms(source_terms=source_terms, source_terms_ids=source_term_ids,
                               target_ontology=EFO_URL, excl_deprecated=True, save_graphs=False,
                               max_mappings=1, min_score=0.0, save_mappings=True, mapper=Mapper.TFIDF,
                               output_file=os.path.join(OUTPUT_FOLDER, f"{dataset_name}_t2t_mappings.csv"),
                               base_iris=("http://www.ebi.ac.uk/efo/", "http://purl.obolibrary.org/obo/MONDO",
                                          "http://purl.obolibrary.org/obo/HP", "http://www.orpha.net/ORDO",
                                          "http://purl.obolibrary.org/obo/DOID"))


def extract_mappings(dataset, dataset_name, entailed_edges_df, source_term_col, source_term_id_col,
                     mapped_term_col, mapped_term_iri_col):
    LOG.info(f"Extracting ontology mappings from {dataset_name}...")
    mappings_list = []
    for _, row in dataset.iterrows():
        mapped_trait_uri = row[mapped_term_iri_col]
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
                if source_term_id_col is None:
                    study_id = onto_utils.generate_uuid()
                else:
                    study_id = row[source_term_id_col]
                mapping = {T2T_INPUT_TERM_ID_COL: study_id,
                           T2T_INPUT_TERM_COL: row[source_term_col],
                           T2T_MAPPED_TERM_COL: row[mapped_term_col],
                           T2T_MAPPED_TERM_IRI_COL: iri,
                           T2T_MAPPED_TERM_CURIE_COL: term_curie,
                           "IS_DISEASE": is_disease,
                           "IS_MEASUREMENT": is_measurement,
                           "IS_PHENOTYPE": is_phenotype,
                           "IS_AMOUNT": is_amount,
                           "IS_PROCESS": is_process}
                mappings_list.append(mapping)
    output_mappings_file = os.path.join(OUTPUT_FOLDER, f"{dataset_name}_mappings.tsv")
    mappings_df = pd.DataFrame(mappings_list)
    mappings_df.to_csv(output_mappings_file, sep="\t", index=False)
    LOG.info(f"...done (saved to {output_mappings_file})")
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


def _fix_curie(value):
    value = value.replace("_", ":")
    if "Orphanet" in value:
        value = value.replace("Orphanet", "ORDO")
    return value


def _get_iri(target_ontology, target_identifier):
    return bioregistry.get_iri(prefix=target_ontology, identifier=target_identifier, priority=['obofoundry', 'default'])


def _get_iri(curie):
    return bioregistry.get_iri(curie, priority=['obofoundry', 'default'])


def compare_mappings(t2t_mappings, benchmark_mappings, edges_df, entailed_edges_df):
    data = []
    queries = t2t_mappings[T2T_INPUT_TERM_ID_COL].unique()

    for query in tqdm(queries):
        t2t_subset = t2t_mappings[t2t_mappings[T2T_INPUT_TERM_ID_COL] == query]
        benchmark_subset = benchmark_mappings[benchmark_mappings[T2T_INPUT_TERM_ID_COL] == query]

        input_trait = t2t_subset[T2T_INPUT_TERM_COL].iloc[0]
        t2t_trait = t2t_subset[T2T_MAPPED_TERM_CURIE_COL].iloc[0]  # there is a single mapping for each input trait
        t2t_trait_label = t2t_subset[T2T_MAPPED_TERM_COL].iloc[0]

        benchmark_traits = benchmark_subset[T2T_MAPPED_TERM_CURIE_COL].unique()
        if len(benchmark_traits) > 0:
            # get asserted EFO parents of both the text2term mapped term and all benchmark mapped terms
            t2t_trait_asserted_parents = edges_df.loc[edges_df['Subject'] == t2t_trait, 'Object'].unique()
            benchmark_traits_asserted_parents = []
            for benchmark_trait in benchmark_traits:
                benchmark_trait_parents = edges_df.loc[edges_df['Subject'] == benchmark_trait, 'Object'].unique()
                benchmark_traits_asserted_parents.extend(benchmark_trait_parents)

            # get entailed parents and children (of the text2term mapped term) in the EFO class hierarchy
            t2t_trait_parents = entailed_edges_df.loc[entailed_edges_df['Subject'] == t2t_trait, 'Object'].unique()
            t2t_trait_children = entailed_edges_df.loc[entailed_edges_df['Object'] == t2t_trait, 'Subject'].unique()

            benchmark_trait = benchmark_traits[0]
            benchmark_trait_label = benchmark_subset[T2T_MAPPED_TERM_COL].iloc[0]

            trait_term = MAPPING_TYPES[benchmark_trait]
            is_disease = trait_term["IS_DISEASE"]
            is_measurement = trait_term["IS_MEASUREMENT"]
            is_phenotype = trait_term["IS_PHENOTYPE"]
            is_amount = trait_term["IS_AMOUNT"]
            is_process = trait_term["IS_PROCESS"]

            if not isinstance(benchmark_trait_label, str):
                try:
                    benchmark_trait_label = benchmark_trait_label[0]
                except TypeError:
                    LOG.error(f"There is no label for benchmark trait {benchmark_trait}.")
                    benchmark_trait_label = ""
            if t2t_trait in benchmark_traits:
                data.append([query, input_trait, t2t_trait, t2t_trait_label, benchmark_trait, benchmark_trait_label,
                             'Same', is_disease, is_measurement, is_phenotype, is_amount, is_process])
            elif any(t in t2t_trait_parents for t in benchmark_traits):
                data.append([query, input_trait, t2t_trait, t2t_trait_label, benchmark_trait, benchmark_trait_label,
                             'More Specific', is_disease, is_measurement, is_phenotype, is_amount, is_process])
            elif any(t in t2t_trait_children for t in benchmark_traits):
                data.append([query, input_trait, t2t_trait, t2t_trait_label, benchmark_trait, benchmark_trait_label,
                             'More General', is_disease, is_measurement, is_phenotype, is_amount, is_process])
            elif any(t in t2t_trait_asserted_parents for t in benchmark_traits_asserted_parents):
                data.append([query, input_trait, t2t_trait, t2t_trait_label, benchmark_trait, benchmark_trait_label,
                             'Sibling', is_disease, is_measurement, is_phenotype, is_amount, is_process])
            else:
                data.append([query, input_trait, t2t_trait, t2t_trait_label, benchmark_trait, benchmark_trait_label,
                             'Unrelated', is_disease, is_measurement, is_phenotype, is_amount, is_process])
        else:
            LOG.debug(f"Could not find benchmark mapping for resource {query}, which was mapped by t2t "
                      f"to '{t2t_trait_label}' ({t2t_trait})")
    return pd.DataFrame(data, columns=[T2T_INPUT_TERM_ID_COL, T2T_INPUT_TERM_COL, "t2t.Mapping", "t2t.MappingLabel",
                                       "Benchmark.Mapping", "Benchmark.MappingLabel", "Classification", "IsDisease",
                                       "IsMeasurement", "IsPhenotype", "IsAmount", "IsProcess"])


def get_gwascatalog_metadata(dataset_name):
    # Load the GWAS Catalog metadata table
    metadata_file = os.path.join("data", "gwascatalog_metadata.tsv")
    LOG.info(f"Loading GWAS Catalog metadata from: {metadata_file}")
    gwascatalog_metadata = pd.read_csv(metadata_file, sep="\t")

    # Filter out the studies/rows that have been mapped to multiple EFO ontology terms
    gwascatalog_metadata = gwascatalog_metadata[
        ~gwascatalog_metadata["MAPPED_TRAIT_CURIE"].astype(str).str.contains(',')]

    benchmark_mappings = extract_mappings(dataset=gwascatalog_metadata,
                                          dataset_name=dataset_name,
                                          entailed_edges_df=efo_entailed_edges_df,
                                          source_term_col="DISEASE.TRAIT",
                                          source_term_id_col="STUDY.ACCESSION",
                                          mapped_term_col="MAPPED_TRAIT",
                                          mapped_term_iri_col="MAPPED_TRAIT_URI")
    return gwascatalog_metadata, benchmark_mappings


def get_efo_ukbb_mappings(dataset_name):
    # Load the UKBB-EFO mappings table
    ukbbefo_table_file = os.path.join("data", "UK_Biobank_master_file.tsv")
    LOG.info(f"Loading UKBB-EFO mappings table from: {ukbbefo_table_file}")
    ukbbefo_table = pd.read_csv(ukbbefo_table_file, sep="\t")

    # Filter out the rows that have been mapped to multiple EFO ontology terms
    ukbbefo_table = ukbbefo_table[~ukbbefo_table["MAPPED_TERM_URI"].astype(str).str.contains(',')]
    ukbbefo_table = ukbbefo_table[~ukbbefo_table["MAPPED_TERM_URI"].astype(str).str.contains(r'\|\|')]

    # Limit table to skos:exactMatch type mappings
    ukbbefo_table = ukbbefo_table[ukbbefo_table["MAPPING_TYPE"] == "Exact"]

    ukbbefo_table["ID"] = [onto_utils.generate_uuid() for _ in range(len(ukbbefo_table))]
    ukbbefo_table['MAPPED_TERM_URI'] = ukbbefo_table['MAPPED_TERM_URI'].apply(_fix_curie)

    benchmark_mappings = extract_mappings(dataset=ukbbefo_table,
                                          dataset_name=dataset_name,
                                          entailed_edges_df=efo_entailed_edges_df,
                                          source_term_col="ZOOMA QUERY",
                                          source_term_id_col="ID",
                                          mapped_term_col="MAPPED_TERM_LABEL",
                                          mapped_term_iri_col="MAPPED_TERM_URI")
    return ukbbefo_table, benchmark_mappings


def get_biomappings(dataset_name, source_ontology="", target_ontology=""):
    # Load the Biomappings table
    biomappings_table_file = os.path.join("data", "biomappings.tsv")
    LOG.info(f"Loading Biomappings table from: {biomappings_table_file}")
    biomappings_table = pd.read_csv(biomappings_table_file, sep="\t")

    # Filter Biomappings table by the source and/or target ontology of mappings
    if source_ontology != "":
        biomappings_table = biomappings_table[biomappings_table['source prefix'] == source_ontology]
    if target_ontology != "":
        biomappings_table = biomappings_table[biomappings_table['target prefix'] == target_ontology]

    # Limit to skos:exactMatch type mappings
    biomappings_table = biomappings_table[biomappings_table['relation'] == 'skos:exactMatch']

    # Add an IRI column to enable IRI-based comparison
    biomappings_table["target IRI"] = biomappings_table.apply(
        lambda row: _get_iri(target_ontology=row["target prefix"], target_identifier=row["target identifier"]), axis=1)

    benchmark_mappings = extract_mappings(dataset=biomappings_table,
                                          dataset_name=dataset_name,
                                          entailed_edges_df=efo_entailed_edges_df,
                                          source_term_col="source name",
                                          source_term_id_col="source identifier",
                                          mapped_term_col="target name",
                                          mapped_term_iri_col="target IRI")
    return biomappings_table, benchmark_mappings


def get_ols_efo_mappings(dataset_name):
    ols_efo_mappings_file = os.path.join("data", "efo.ols.sssom.tsv")
    LOG.info(f"Loading OLS mappings table from: {ols_efo_mappings_file}")
    ols_efo_mappings = pd.read_csv(ols_efo_mappings_file, sep="\t", skiprows=112)

    # Limit table to rows whose mappings have a label to work with
    ols_efo_mappings = ols_efo_mappings.dropna(subset=['object_label'])

    ols_efo_mappings["efo_IRI"] = ols_efo_mappings['subject_id'].apply(_get_iri)
    # ols_efo_mappings["target IRI"] = ols_efo_mappings.apply(lambda row: _get_iri(row["subject_id"]), axis=1)

    benchmark_mappings = extract_mappings(dataset=ols_efo_mappings,
                                          dataset_name=dataset_name,
                                          entailed_edges_df=efo_entailed_edges_df,
                                          source_term_col="object_label",
                                          source_term_id_col="object_id",
                                          mapped_term_col="subject_label",
                                          mapped_term_iri_col="efo_IRI")
    return ols_efo_mappings, benchmark_mappings


if __name__ == '__main__':
    # Benchmark mapping sets
    gc = "GWASCatalog"
    ub = "UKBB-EFO"
    bm = "Biomappings"
    oe = "OLS-EFO"

    # Create output directory if it does not exist
    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)

    # Load the EFO ontology table containing all asserted SubClassOf relationships between terms in EFO
    efo_edges_df = pd.read_csv(os.path.join("data", "efo_edges.tsv"), sep="\t")

    # Load the EFO ontology table containing all entailed SubClassOf relationships between terms in EFO
    efo_entailed_edges_df = pd.read_csv(os.path.join("data", "efo_entailed_edges.tsv"), sep="\t")

    # Compare text2term mappings to the UK Biobank-EFO mapping set
    efo_ukbb, efo_ukbb_mappings = get_efo_ukbb_mappings(ub)
    compare(benchmark_dataset=efo_ukbb, benchmark_mappings=efo_ukbb_mappings, benchmark_dataset_name=ub,
            source_term_col="ZOOMA QUERY", source_term_id_col="ID")

    # Compare text2term mappings to the Biomappings set of EFO mappings
    biomappings, biomappings_efo_mappings = get_biomappings(bm, target_ontology="efo")
    compare(benchmark_dataset=biomappings, benchmark_mappings=biomappings_efo_mappings, benchmark_dataset_name=bm,
            source_term_col="source name", source_term_id_col="source identifier")

    # Compare text2term mappings to the GWAS Catalog mapping set
    gwascatalog, gwascatalog_mappings = get_gwascatalog_metadata(gc)
    compare(benchmark_dataset=gwascatalog, benchmark_mappings=gwascatalog_mappings, benchmark_dataset_name=gc,
            source_term_col="DISEASE.TRAIT", source_term_id_col="STUDY.ACCESSION")

    # Compare text2term mappings to the OLS mappings extracted from EFO
    ols_efo, ols_mappings = get_ols_efo_mappings(oe)
    compare(benchmark_dataset=ols_efo, benchmark_mappings=ols_mappings, benchmark_dataset_name=oe,
            source_term_col="object_label", source_term_id_col="object_id")
