from typing import List, Set, Dict, Any
import numpy


# Main ontology to reduced dag functions
# --------------------------------------------------------------------------------------------------
def ontology_to_weighted_dag(concepts, abundances, root, ontology_dag, show_lvs):
    classified_concepts = classify_concepts(concepts, ontology_dag)
    concepts_all_classes = get_all_classes(classified_concepts, ontology_dag, root)
    abundances_dict = get_abundance_dict(abundances, concepts)
    calculated_weights = calculate_weights(concepts_all_classes, abundances_dict, show_lvs)
    return calculated_weights


def reduce_d_ontology(complete_dictionary: Dict[str, Any],
                      classes_abundance: Dict[str, float]) -> Dict[str, Any]:
    """ Extract the sub-graph of the d_classes_ontology dictionary conserving only nodes implicated
    with the concepts studied.

    Parameters
    ----------
    complete_dictionary: Dict[str, Any]
        Dictionary of the ontology complete graph
    classes_abundance: Dict[str, float]
        Dictionary of abundances (keys are all nodes implicated to be conserved)

    Returns
    -------
    Dict[str, Any]
        Dictionary of the ontology sub-graph conserving only nodes implicated with the concepts
        studied.
    """
    if complete_dictionary is not None:
        reduced_dictionary = dict()
        for k, v in complete_dictionary.items():
            if k in classes_abundance:
                reduced_dictionary[k] = v
        return reduced_dictionary

# ==================================================================================================
# REDUCE DAG FUNCTIONS
# ==================================================================================================

def classify_concepts(concepts: List[str], ontology_dag: Dict[str, List[str]]) \
        -> Dict[str, List[str]]:
    """ Extract concepts able to be classified in the ontology.

    Parameters
    ----------
    concepts: List[str]
        List of concepts to classify
    ontology_dag: Dict[str, List[str]]
        Dictionary of the classes ontology associating for each concept its +1 parent classes.

    Returns
    -------
    Dict[str, List[str]]
        Dictionary associating to each classified concept all its +1 parent classes.
    """
    classified_concepts = dict()
    print(f'{len(concepts)} concepts to classify')
    for cpt in concepts:
        try:
            classified_concepts[cpt] = ontology_dag[cpt]
            classified = True
        except KeyError:
            classified = False
        if not classified:
            print(f'{cpt} not classified.')
    print(f'{len(classified_concepts)}/{len(concepts)} concepts classified')
    return classified_concepts


# Recursive class extraction function
# --------------------------------------------------------------------------------------------------
def get_all_classes(obj_classes: Dict[str, List[str]], d_classes_ontology: Dict[str, List[str]],
                    root_item: str) -> Dict[str, Set[str]]:
    """ Extract all parent classes for each metabolite.

    Parameters
    ----------
    obj_classes: Dict[str, List[str]] (Dict[metabolite, List[class]])
        Dictionary associating for each object the list of +1 parent classes it belongs to.
    d_classes_ontology: Dict[str, List[str]]
        Dictionary of the classes ontology associating for each class its +1 parent classes.
    root_item: str
        Name of the root item of the ontology.

    Returns
    -------
    Dict[str, Set[str]] (Dict[metabolite, Set[class]])
        Dictionary associating for each metabolite the list of all parent classes it belongs to.
    """
    all_classes_met = dict()
    for met, classes in obj_classes.items():
        all_classes = set(classes)
        for c in classes:
            if c != root_item:
                m_classes = get_parents(c, set(d_classes_ontology[c]), d_classes_ontology, root_item)
                all_classes = all_classes.union(m_classes)
        all_classes_met[met] = all_classes
    return all_classes_met


def get_parents(child: str, parent_set: Set[str], d_classes_ontology: Dict[str, List[str]],
                root_item) -> Set[str]:
    """ Get recursively from a child class, all its parents classes found in ontology.

    Parameters
    ----------
    child: str
        Child class
    parent_set: Set[str]
        Set of all parents from previous classes
    d_classes_ontology: Dict[str, List[str]]
        Dictionary of the classes ontology of MetaCyc associating for each class its parent classes.
    root_item: str
        Name of the root item of the ontology

    Returns
    -------
    Set[str]
        Set of the union of the set  of child parent classes and the set of all previous parents.
    """
    parents = d_classes_ontology[child]
    for p in parents:
        parent_set.add(p)
    for p in parents:
        if p != root_item:
            parent_set = get_parents(p, parent_set, d_classes_ontology, root_item)
    return parent_set


# ==================================================================================================
# WEIGHTS CALCULATION
# ==================================================================================================

def get_abundance_dict(abundances: List[float] or None, concepts: List[str])\
        -> Dict[str, float]:
    """ Generate abundances dictionary.

    Parameters
    ----------
    abundances: List[float] (size N) or None
        List of concepts abundances (or None if no abundances associated --> will associate
        an abundance of 1 for each concept)
    concepts: List[str] (size N)
        List of concepts ID.

    Returns
    -------
    Dict[str, float]
        Dictionary associating to each concept its abundance.
    """
    if abundances is None:
        abundances = len(concepts) * [1]
    if len(concepts) == len(abundances):
        abundances_dict = {}
        for i in range(len(concepts)):
            abundances_dict[concepts[i]] = abundances[i]
    else:
        raise AttributeError(f'Length of concepts IDs list must be equal to '
                             f'its abundances list length : {len(concepts)} '
                             f'!= {len(abundances)}')
    return abundances_dict


def calculate_weights(all_classes: Dict[str, Set[str]], abundances_dict: Dict[str, float],
                      show_leaves: bool) -> Dict[str, float]:
    """ Indicate for each class the number of base object found belonging to the class

    Parameters
    ----------
    all_classes: Dict[str, Set[str]] (Dict[metabolite, Set[class]])
        Dictionary associating for each concept the list of all parent classes it belongs to.
    abundances_dict: Dict[str, float]
        Dictionary associating for each concept, its abundance value
    show_leaves: bool
        True to show input metabolic objets at sunburst leaves

    Returns
    -------
    Dict[str, float]
        Dictionary associating for each class the weight of concepts found belonging to the class.
    """
    classes_abondance = dict()
    for met, classes in all_classes.items():
        if show_leaves:
            if met not in classes_abondance.keys():
                classes_abondance[met] = abundances_dict[met]
            else:
                classes_abondance[met] += abundances_dict[met]
        for c in classes:
            if c not in classes_abondance.keys():
                classes_abondance[c] = abundances_dict[met]
            else:
                classes_abondance[c] += abundances_dict[met]
    return dict(reversed(sorted(classes_abondance.items(), key=lambda item: item[1])))


def get_classes_scores(all_classes, scores_dict, root):
    if scores_dict is not None:
        classes_scores = dict()
        for met, classes in all_classes.items():
            if met in scores_dict.keys():
                classes_scores[met] = scores_dict[met]
            else:
                classes_scores[met] = numpy.nan
        classes_scores[root] = numpy.nan
        return classes_scores
