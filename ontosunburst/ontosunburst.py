import os
import json
from typing import List, Dict, Set
from time import time
import plotly.graph_objects as go

from ontosunburst.onto2dag import ontology_to_weighted_dag, get_classes_scores, reduce_d_ontology


from ontosunburst.dag2tree import TreeData, get_name, BINOMIAL_TEST, HYPERGEO_TEST, ROOT_CUT, \
    ROOT_TOTAL_CUT, ROOT_UNCUT, PATH_UNCUT, PATH_BOUND, PATH_DEEPER, PATH_HIGHER
from ontosunburst.tree2sunburst import generate_sunburst_fig, TOPOLOGY_A, ENRICHMENT_A

# ==================================================================================================
#                                           CONSTANTS
# ==================================================================================================

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
DEFAULT_PATH = os.path.join(CURRENT_DIR, 'Inputs')
CLASSES_SUFFIX = 'classes.json'
LABELS_SUFFIX = 'labels.json'

METACYC = 'metacyc'
EC = 'ec'
CHEBI = 'chebi'
CHEBI_R = 'chebi_r'
GO_CC = 'go_cc'
GO_MF = 'go_mf'
GO_BP = 'go_bp'
GO = 'go'
KEGG = 'kegg'

ROOTS = {METACYC: 'FRAMES',
         CHEBI: 'chebi',
         CHEBI_R: 'CHEBI:50906',
         EC: 'Enzyme',
         GO_CC: 'GO:0005575',
         GO_BP: 'GO:0008150',
         GO_MF: 'GO:0003674',
         GO: 'GO',
         KEGG: 'kegg'}


# ==================================================================================================
#                                            WORKFLOW
# ==================================================================================================

def ontosunburst(interest_set: List[str],
                 ontology: str = None,
                 abundances: List[float] = None,
                 reference_set: List[str] = None,
                 ref_abundances: List[float] = None,
                 analysis: str = TOPOLOGY_A,
                 output: str = 'sunburst',
                 scores: Dict[str, float] = None,
                 write_output: bool = True,
                 ontology_dag_input: str or Dict[str, str] = None,
                 input_root: str = None,
                 id_to_label_input: str or Dict[str, str] = None,
                 labels: bool = True,
                 test: str = BINOMIAL_TEST,
                 root_cut: str = ROOT_CUT,
                 path_cut: str = PATH_UNCUT,
                 ref_base: bool = False,
                 show_leaves: bool = False,
                 **kwargs) -> go.Figure:
    """ Main function to be called generating the sunburst figure

    Parameters
    ----------
    interest_set: List[str]
        Interest list of concepts IDs to classify.
    ontology: str (optional, default=None, values in ['metacyc', 'ec', 'chebi', 'chebi_r', 'kegg',
                                                      'go_cc', 'go_bp', 'go_mf', 'go', None])
        Ontology name to use.
    abundances: List[str] (optional, default=None)
        Abundance values associated to interest_set list parameter
    reference_set: List[str] (optional, default=None)
        Reference list of concepts IDs.
    ref_abundances: List[str] (optional, default=None)
        Abundance values associated to reference_set list parameter
    analysis: str (optional, default='topology', values in ['topology', 'enrichment'])
        Analysis mode : topology or enrichment.
    output: str (optional, default='sunburst')
        Path of the output to save figure, if None, outputs will be sunburst.html and sunburst.tsv
        files
    scores: Dict[str, float] (optional, default=None)
        Dictionary associating for each ontology ID, its precalculated enrichment score. If None
        enrichment will be calculated.
    write_output: bool (optional, default=True)
        True to write the html figure and tsv class files, False to only return plotly sunburst
        figure.
    ontology_dag_input: str or Dict[str, str] (optional, default=None)
        Ontology DAG dictionary or json file. Use if tailored ontology or alternative
        (modified, updated, ...) default ontology DAG.
    input_root: str (optional, default=None)
        Root item of the ontology  (to precise if tailored ontology).
    id_to_label_input: str or Dict[str, str] (optional, default=None)
        Path to ID-LABELS association json file or ID-LABELS association dictionary.
        If None default files will be used. Use if tailored ontology or alternative
        (modified, updated, ...) default ontology.
    labels: bool (optional, default=True)
        True to show labels as sunburst sectors labels, False to show ID as sunburst sectors labels.
    test: str (optional, default='binomial', values in ['binomial', 'hypergeometric'])
        Type of test if analysis=enrichment, binomial or hypergeometric test.
    root_cut: str (optional, default='cut', values in ['uncut', 'cut', 'total'])
        mode for root cutting (uncut, cut or total)
    path_cut: str (optional, default='uncut', values in ['uncut', 'deeper', 'higher', 'bound'])
        mode for nested path cutting (uncut, deeper, higher or bound)
    ref_base: bool (optional, default=False)
        True to have the base classes representation of the reference set in the figure.
    show_leaves: bool (optional, default=False)
        True to show input metabolic objets at sunburst leaves
    **kwargs

    Returns
    -------
    go.Figure
        Plotly graph_objects figure of the sunburst
    """
    start_time = time()
    # LOAD ID TO LABELS DICTIONARY -----------------------------------------------------------------
    id_to_label = get_id_to_label_dict(id_to_label_input, labels, ontology)
    # LOAD ONTOLOGY DAG DICTIONARY -----------------------------------------------------------------
    ontology_dag = get_ontology_dag_dict(ontology, ontology_dag_input)
    # GET ROOT -------------------------------------------------------------------------------------
    root = get_ontology_root(ontology, input_root)
    # WORKFLOW -------------------------------------------------------------------------------------
    fig = _global_analysis(analysis=analysis,
                           interest_concepts=interest_set, abundances=abundances,
                           scores=scores,
                           reference_concepts=reference_set, ref_abundances=ref_abundances,
                           ontology_dag=ontology_dag,
                           output=output, write_output=write_output, id_to_label=id_to_label,
                           test=test, root=root, root_cut=root_cut, path_cut=path_cut,
                           ref_base=ref_base, show_leaves=show_leaves, **kwargs)
    end_time = time()
    print(f'Execution time : {end_time - start_time} seconds')
    return fig


def _global_analysis(analysis, interest_concepts, abundances, scores, reference_concepts,
                     ref_abundances, ontology_dag, output, write_output, id_to_label,
                     test, root, root_cut, path_cut, ref_base, show_leaves, **kwargs):
    """

    Parameters
    ----------
    analysis
    interest_concepts
    abundances
    reference_concepts
    ref_abundances
    ontology_dag
    output
    write_output
    id_to_label
    test
    root
    root_cut
    path_cut
    ref_base
    show_leaves
    kwargs

    Returns
    -------

    """
    # ONTOLOGY TO WEIGHTED DAG
    # =============================================================================================
    # Calculate all concepts weights --------------------------------------------------------------
    calculated_weights = ontology_to_weighted_dag(concepts=interest_concepts, abundances=abundances,
                                                  root=root, ontology_dag=ontology_dag,
                                                  show_lvs=show_leaves)

    if reference_concepts is not None:
        ref_set = True
        ref_calculated_weights = ontology_to_weighted_dag(concepts=reference_concepts,
                                                          abundances=ref_abundances, root=root,
                                                          ontology_dag=ontology_dag,
                                                          show_lvs=show_leaves)
    else:
        ref_set = False
        ref_calculated_weights = calculated_weights

    # Scores
    classes_scores = get_classes_scores(calculated_weights, scores, root)

    # Reduce ontology (get DAG subgraph) ----------------------------------------------------------
    if ref_base:
        ontology_dag = reduce_d_ontology(ontology_dag, ref_calculated_weights)
        id_to_label = reduce_d_ontology(id_to_label, ref_calculated_weights)
    else:
        ontology_dag = reduce_d_ontology(ontology_dag, calculated_weights)
        id_to_label = reduce_d_ontology(id_to_label, calculated_weights)

    # WRITE CONCEPTS CLASSES IN TSV OUTPUT FILE ---------------------------------------------------
    # if write_output:
    #     write_concepts_classes(ontology, concepts_all_classes, output, id_to_label)

    # DAG TO TREE
    # =============================================================================================
    tree_data = TreeData()
    tree_data.dag_to_tree(set_abundance=calculated_weights, ref_abundance=ref_calculated_weights,
                          parent_dict=ontology_dag, root_item=root, names=id_to_label,
                          ref_base=ref_base)

    tree_data.calculate_proportions(ref_base)
    significant = None
    if analysis == ENRICHMENT_A:
        significant = tree_data.make_enrichment_analysis(test, classes_scores)
    tree_data.cut_root(root_cut)
    tree_data.cut_nested_path(path_cut, ref_base)

    # TREE TO SUNBURST
    # =============================================================================================
    return generate_sunburst_fig(data=tree_data, output=output, analysis=analysis, test=test,
                                 significant=significant, ref_set=ref_set,
                                 write_fig=write_output, **kwargs)


# ==================================================================================================
#                                             FUNCTIONS
# ==================================================================================================

def get_file(ontology, suffix):
    for file in os.listdir(DEFAULT_PATH):
        if file.startswith(ontology + '__') and file.endswith('__' + suffix):
            return os.path.join(DEFAULT_PATH, file)


def aggregate_go_ontologies(suffix):
    go_aggregated = dict()
    for sub_go_ontology in [GO_BP, GO_CC, GO_MF]:
        dict_sub_onto_input = get_file(sub_go_ontology, suffix)
        with open(dict_sub_onto_input, 'r') as f:
            dict_sub_onto = json.load(f)
        go_aggregated.update(dict_sub_onto)
    if suffix == CLASSES_SUFFIX:
        for sub_go_ontology in [GO_BP, GO_CC, GO_MF]:
            go_aggregated[ROOTS[sub_go_ontology]] = [ROOTS[GO]]
    return go_aggregated


def get_id_to_label_dict(id_to_label_input, labels, ontology):
    # Returns ID_to_Labels only if labels is True
    if labels:
        # Case default ontology AND use of default labels file
        if ontology is not None and id_to_label_input is None:
            if ontology == GO:
                return aggregate_go_ontologies(LABELS_SUFFIX)
            id_to_label_input = get_file(ontology, LABELS_SUFFIX)
        # Case id_to_label_input parameter filled
        if id_to_label_input is not None:
            # Case id_to_label_input parameter is a file path (str)
            if type(id_to_label_input) == str:
                with open(id_to_label_input, 'r') as f:
                    id_to_label = json.load(f)
                    return id_to_label
            # Case id_to_label_input parameter is a dictionary (dict)
            elif type(id_to_label_input) == dict:
                return id_to_label_input
            # Case id_to_label_input parameter is not a dictionary (dict), neither a file
            # path (str) : raises an error
            else:
                raise ValueError('id_to_label_input parameter must be a json file path (str) or a '
                                 'dictionary')


def get_ontology_dag_dict(ontology, ontology_dag_input):
    # Case ontology_dag_input parameter not filled (default : None)
    if ontology_dag_input is None:
        # Case no default ontology : raises an error
        if ontology is None:
            raise ValueError('If no default ontology, must fill ontology_dag_input parameter')
        # Case default ontology : get default ontology file path
        else:
            if ontology == GO:
                return aggregate_go_ontologies(CLASSES_SUFFIX)
            ontology_dag_input = get_file(ontology, CLASSES_SUFFIX)
    # Case ontology_dag_input parameter is a file path (str)
    if type(ontology_dag_input) == str:
        with open(ontology_dag_input, 'r') as f:
            ontology_dag = json.load(f)
            return ontology_dag
    # Case ontology_dag_input parameter is a dictionary (dict)
    elif type(ontology_dag_input) == dict:
        return ontology_dag_input
    # Case ontology_dag_input parameter is not a dictionary (dict), neither a file path (str) :
    # raises an error
    else:
        raise ValueError('ontology_dag_input parameter must be a json file path (str) or a '
                         'dictionary')


def get_ontology_root(ontology, input_root):
    if ontology is not None:
        return ROOTS[ontology]
    elif input_root is None:
        raise ValueError('If no default ontology, must fill input_root parameter')
    else:
        return input_root


def write_concepts_classes(ontology: str, all_classes: Dict[str, Set[str]], output: str,
                           id_to_label: Dict[str, str]):
    """ Writes, for each input class, all its ancestors in a .tsv file.

    Parameters
    ----------
    ontology
    all_classes
    output
    id_to_label
    """
    if ontology is None:
        ontology = ''
    links_dict = {METACYC: 'https://metacyc.org/compound?orgid=META&id=',
                  CHEBI: 'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:',
                  CHEBI_R: 'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:',
                  EC: 'https://enzyme.expasy.org/EC/',
                  KEGG: 'https://www.genome.jp/entry/',
                  GO_MF: 'https://amigo.geneontology.org/amigo/term/',
                  GO_CC: 'https://amigo.geneontology.org/amigo/term/',
                  GO_BP: 'https://amigo.geneontology.org/amigo/term/',
                  GO: 'https://amigo.geneontology.org/amigo/term/',
                  '': ''}
    with open(f'{output}.tsv', 'w') as f:
        f.write('\t'.join(['ID', 'Label', 'Classes ID', 'Classes Label', 'Link']) + '\n')
        for met_id, classes_id, in all_classes.items():
            link = links_dict[ontology] + met_id
            met_lab = get_name(met_id, id_to_label)
            classes_lab = [get_name(cl, id_to_label) for cl in classes_id]
            f.write('\t'.join([met_id, met_lab, ', '.join(classes_id), ', '.join(classes_lab),
                               link]) + '\n')