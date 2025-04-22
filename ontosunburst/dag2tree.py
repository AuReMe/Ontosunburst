from typing import List, Dict, Set
import numpy as np
from numpy import nan
import scipy.stats as stats

# ==================================================================================================
# CONSTANTS
# ==================================================================================================

# Comparison tests
BINOMIAL_TEST = 'binomial'
HYPERGEO_TEST = 'hypergeometric'

# Analysis method
TOPOLOGY_A = 'topology'
ENRICHMENT_A = 'enrichment'

MAX_RELATIVE_NB = 1000000

# Keys
# ----
IDS = 'ID'
ONTO_ID = 'Onto ID'
PARENT = 'Parent'
LABEL = 'Label'
WEIGHT = 'Weight'
REF_WEIGHT = 'Reference weight'
PROP = 'Proportion'
REF_PROP = 'Reference proportion'
RELAT_PROP = 'Relative proportion'
PVAL = 'Pvalue'

# Root cut
ROOT_CUT = 'cut'
ROOT_TOTAL_CUT = 'total'
ROOT_UNCUT = 'uncut'

# Path cut
PATH_UNCUT = 'uncut'
PATH_DEEPER = 'deeper'
PATH_HIGHER = 'higher'
PATH_BOUND = 'bound'


# ==================================================================================================
# CLASS
# ==================================================================================================
class TreeData:
    C_ID = 0
    """
    TreeData class: stores figure parameters values.

    Attributes
    ----------
    self.ids: List[str]
        Unique ids of all sectors (id=path in the tree)
    self.onto_ids: List[str]
        Ontology ID of all sectors
    self.labels: List[str]
        Sectors labels
    self.parents: List[str]
        Sectors parents
    self.count: List[float]
        Sectors interest set count (0<count)
    self.ref_count: List[float]
        Sectors reference set count (0<count)
    self.prop: List[float]
        Sectors interest set proportion (0<prop<1)
    self.ref_prop: List[float]
        Sectors reference set proportion (0<prop<1)
    self.relative_prop: List[int]
        Sectors relative proportion (0<r_prop<1000000)
    self.p_val: List[float]
        Sectors p-value if enrichment analysis
    self.len: int
        Number of sectors
    """

    def __init__(self):
        self.ids = list()
        self.onto_ids = list()
        self.labels = list()
        self.parents = list()
        self.count = list()
        self.ref_count = list()
        self.prop = list()
        self.ref_prop = list()
        self.relative_prop = list()
        self.p_val = list()
        self.len = 0

    def __str__(self):
        string = ''
        data = self.get_data_dict()
        for k, v in data.items():
            string += f'{k}\n{"-" * len(k)}\n{v}\n'
        return string

    def get_data_dict(self):
        return {IDS: self.ids, ONTO_ID: self.onto_ids, LABEL: self.labels, PARENT: self.parents,
                WEIGHT: self.count, REF_WEIGHT: self.ref_count, PROP: self.prop,
                REF_PROP: self.ref_prop, RELAT_PROP: self.relative_prop, PVAL: self.p_val}

    def dag_to_tree(self, set_abundance: Dict[str, float], ref_abundance: Dict[str, float],
                    parent_dict: Dict[str, List[str]], root_item: str,
                    names: Dict[str, str] = None, ref_base: bool = True):
        """ Fill TreeData list attributes (self.ids, self.onto_ids, self.labels, self.parents,
        self.count, self.ref_count)

        Parameters
        ----------
        set_abundance: Dict[str, float]
            Dictionary associating for each class the number of objects found belonging to the class
            in the interest set
        ref_abundance: Dict[str, float]
            Dictionary associating for each class the number of objects found belonging to the class
            in the reference set
        parent_dict: Dict[str, List[str]]
            Dictionary associating for each class, its parents classes
        root_item: str
            Name of the root item of the ontology
        names: Dict[str, str]
            Dictionary associating for some or each ontology IDs, its label
        ref_base: bool
            True to have the reference as base, False otherwise
        """
        children_dict = get_children_dict(parent_dict)
        self.dag_traversal_rec(root_item, children_dict, names, ref_abundance, set_abundance,
                               ref_base, '')

    def dag_traversal_rec(self, c_onto_id: str, children_dict: Dict[str, List[str]],
                          names: Dict[str, str], ref_abundance: Dict[str, float],
                          set_abundance: Dict[str, float], ref_base: bool, p_id: str):
        """ Fill parameters recursively from the root. Perform a traversing of the DAG and create a
        vertex of a tree for each visited node (even if already visited, in this case vertex are
        duplicated with the same label but a different ID)

        Parameters
        ----------
        c_onto_id: str
            Ontology ID of the current concept visited
        children_dict: Dict[str, List[str]]
            Dictionary associating for each concept, the list of its -1 children concepts
        names: Dict[str, str]
            Dictionary associating for some or each ontology IDs, its label
        ref_abundance: Dict[str, float]
            Dictionary associating for each class the number of objects found belonging to the class
            in the reference set
        set_abundance: Dict[str, float]
            Dictionary associating for each class the number of objects found belonging to the class
            in the interest set
        ref_base: bool
            True to have the reference as base, False otherwise
        p_id: str
            ID (not ontology ID) of the parent concept
        """
        if (ref_base and c_onto_id in ref_abundance) or \
                (not ref_base and c_onto_id in set_abundance):
            self.C_ID += 1
            c_id = str(self.C_ID)
            c_label = get_name(c_onto_id, names)
            c_ref_abundance = ref_abundance[c_onto_id]
            c_abundance = get_set2_abundance(set_abundance, c_onto_id)
            self.add_value(m_id=c_id, onto_id=c_onto_id, label=c_label, count=c_abundance,
                           ref_count=c_ref_abundance, parent=p_id)
            for child in children_dict[c_onto_id]:
                self.dag_traversal_rec(child, children_dict, names, ref_abundance, set_abundance,
                                       ref_base, c_id)

    def add_value(self, m_id: str, onto_id: str, label: str, count: float, ref_count: float,
                  parent: str):
        """ Fill the data attributes for an object class.

        Parameters
        ----------
        m_id: str
            ID unique of the object class to add
        onto_id: str
            ID in the ontology
        label: str
            Label (name) of the object class to add
        count: int
            Abundance value of the object class to add
        ref_count: int
            Reference abundance value of the object class to add
        parent: str
            Parent object class of the object class to add
        """
        if m_id in self.ids:
            raise ValueError(f'{m_id} already in data IDs, all IDs must be unique.')
        self.ids.append(m_id)
        self.onto_ids.append(onto_id)
        self.labels.append(label)
        self.parents.append(parent)
        self.count.append(count)
        self.ref_count.append(ref_count)
        self.prop.append(nan)
        self.ref_prop.append(nan)
        self.relative_prop.append(nan)
        self.p_val.append(nan)
        self.len += 1

    def calculate_proportions(self, ref_base: bool):
        """ Calculate TreeData proportion list attributes (self.prop, self.ref_prop,
        self.relative_prop). If total add relative proportion to +1 parent for branch value.

        Parameters
        ----------
        ref_base: bool
            True if reference base representation
        """
        # Get total proportion
        max_abondance = int(np.nanmax(self.count))
        self.prop = [x / max_abondance for x in self.count]
        # Get reference proportion
        max_ref_abondance = np.max(self.ref_count)
        self.ref_prop = [x / max_ref_abondance for x in self.ref_count]
        # Get proportion relative to +1 parent proportion for total branch value
        self.relative_prop = [x for x in self.prop]
        p = ''
        self.__get_relative_prop(p, ref_base)
        # IDK WHY IT WORKS ???
        missed = [self.ids[i] for i in range(self.len) if self.relative_prop[i] < 1]
        if missed:
            parents = {self.parents[self.ids.index(m)] for m in missed}
            for p in parents:
                self.__get_relative_prop(p, ref_base)

    def __get_relative_prop(self, p_id: str, ref_base: bool):
        """ Get recursively relative proportion of a parent children to itself. Set it to class
        self.relative_prop attribute.

        Parameters
        ----------
        p_id: str
            ID of the parent
        ref_base: bool
            True if reference base representation
        """
        if ref_base:
            base_count = self.ref_count
        else:
            base_count = self.count
        if p_id == '':
            prop_p = MAX_RELATIVE_NB
            count_p = max(base_count)
        else:
            prop_p = self.relative_prop[self.ids.index(p_id)]
            count_p = base_count[self.ids.index(p_id)]
        index_p = [i for i, v in enumerate(self.parents) if v == p_id]
        p_children = [self.ids[i] for i in index_p]
        count_p_children = [base_count[i] for i in index_p]
        if np.nansum(count_p_children) > count_p:
            total = np.nansum(count_p_children)
        else:
            total = count_p
        for i, c in enumerate(p_children):
            if not ref_base and np.isnan(self.prop[self.ids.index(c)]):
                prop_c = 0
            else:
                prop_c = int((count_p_children[i] / total) * prop_p)
            self.relative_prop[self.ids.index(c)] = prop_c
        for c in p_children:
            if c in self.parents:
                self.__get_relative_prop(c, ref_base)

    def make_enrichment_analysis(self, test: str, scores: Dict[str, float] = None) \
            -> Dict[str, float]:
        """ Performs statistical tests for enrichment analysis.

        Parameters
        ----------
        test: str
            Type of test : binomial or hypergeometric
        scores: Dict[str, float]
            Dictionary associating for each ontology ID, its enrichment score. If None enrichment
            will be calculated.

        Returns
        -------
        Dict[str, float]
            Dictionary of significant metabolic object label associated with their p-value
        """
        nb_classes = len(set([self.labels[i] for i in range(self.len)
                              if not np.isnan(self.count[i])]))
        significant_representation = dict()
        if scores is not None:
            for i in range(self.len):
                p_val = scores[self.onto_ids[i]]
                self.p_val[i] = -np.log10(p_val)
                if p_val < 0.05 / nb_classes:  # Keep significant p-values : Bonferroni
                    significant_representation[self.onto_ids[i]] = p_val
        else:
            m = np.max(self.ref_count)  # M = ref set total item number
            n = int(np.nanmax(self.count))  # N = interest set total item number
            for i in range(self.len):
                if type(self.count[i]) == int:  # If count not nan (= if concept in interest set)
                    # Binomial Test
                    if test == BINOMIAL_TEST:
                        p_val = stats.binomtest(self.count[i], n, self.ref_count[i] / m,
                                                alternative='two-sided').pvalue
                    # Hypergeometric Test
                    elif test == HYPERGEO_TEST:
                        p_val_upper = stats.hypergeom.sf(self.count[i] - 1, m, self.ref_count[i], n)
                        p_val_lower = stats.hypergeom.cdf(self.count[i], m, self.ref_count[i], n)
                        p_val = 2 * min(p_val_lower, p_val_upper)  # bilateral
                    else:
                        raise ValueError(
                            f'test parameter must be in : {[BINOMIAL_TEST, HYPERGEO_TEST]}')
                    if ((self.count[i] / n) - (self.ref_count[i] / m)) > 0:  # If over-represented :
                        self.p_val[i] = -np.log10(p_val)  # Positive log10(p-value)
                    else:  # If under-represented :
                        self.p_val[i] = np.log10(p_val)  # Negative log10(p-value)
                    if p_val < 0.05 / nb_classes:  # Keep significant p-values : Bonferroni
                        significant_representation[self.onto_ids[i]] = p_val
        significant_representation = dict(
            sorted(significant_representation.items(), key=lambda item: item[1]))
        return significant_representation

    def cut_root(self, mode: str):
        """ Filter data to cut (or not) the root to remove not necessary 100% represented classes.

        Parameters
        ----------
        mode: str
            Mode of root cutting
            - uncut: doesn't cut and keep all nodes from ontology root
            - cut: keep only the lowest level 100% shared node
            - total: remove all 100% shared nodes (produces a pie at center)
        """
        if mode not in {ROOT_UNCUT, ROOT_CUT, ROOT_TOTAL_CUT}:
            raise ValueError(f'Root cutting mode {mode} unknown, '
                             f'must be in {[ROOT_UNCUT, ROOT_CUT, ROOT_TOTAL_CUT]}')
        if mode == ROOT_CUT or mode == ROOT_TOTAL_CUT:
            roots_ind = [i for i in range(self.len) if self.relative_prop[i] == MAX_RELATIVE_NB]
            roots = [self.ids[i] for i in roots_ind]
            roots_lab = [self.labels[i] if self.labels[i] not in self.ids
                         else self.labels[i] + '_' for i in roots_ind]
            lab = {roots[i]: roots_lab[i] for i in range(len(roots))}
            self.delete_value(roots_ind)
            if mode == ROOT_CUT:
                self.parents = [lab[x] if x in roots else x for x in self.parents]
            if mode == ROOT_TOTAL_CUT:
                self.parents = ['' if x in roots else x for x in self.parents]

    def cut_nested_path(self, mode: str, ref_base: bool):
        """ Cut nested path in the tree graph (path of nested sectors sharing the same value)

        Parameters
        ----------
        mode: str
            Mode of path cutting
            - uncut: doesn't cut and keep all sectors
            - deeper: cut nested path and only conserve the deepest sector in the tree
            - higher: cut nested path and only conserve the highest sector in the tree
            - bound: cut nested path and only conserve the highest AND the deepest sectors in the
            tree
        ref_base: bool
            True if reference base representation
        """
        if ref_base:
            count = self.ref_count
        else:
            count = self.count
        if mode != PATH_UNCUT:
            nested_paths = []
            for p_i in range(self.len):
                p = self.ids[p_i]
                p_children = [self.ids[i] for i in range(self.len) if self.parents[i] == p]
                if len(p_children) == 1:
                    p_p = self.parents[p_i]
                    p_p_children = [self.ids[i] for i in range(self.len) if self.parents[i] == p_p]
                    if len(p_p_children) != 1:
                        p_count = count[p_i]
                        c_i = self.ids.index(p_children[0])
                        c_count = count[c_i]
                        if p_count == c_count:
                            nested_paths.append(self.get_full_nested_path(c_i, [p_i], count))
            self.delete_nested_path(mode, nested_paths)

    def get_full_nested_path(self, p_i: int, n_path: List[int], count: List[float]):
        """ Get all index of a nested path sector from its parent sector index.

        Parameters
        ----------
        p_i: int
            Parent sector index of the nested path
        n_path: List[int]
            List of sector indexes of the nested path
        count: List[float]
            List of all sectors count value.

        Returns
        -------
        List[int]
            List of sector indexes of the nested path
        """
        n_path.append(p_i)
        p = self.ids[p_i]
        p_children = [self.ids[i] for i in range(self.len) if self.parents[i] == p]
        if len(p_children) == 1:
            p_count = count[p_i]
            c_i = self.ids.index(p_children[0])
            c_count = count[c_i]
            if p_count == c_count:
                n_path = self.get_full_nested_path(c_i, n_path, count)
        return n_path

    def delete_nested_path(self, mode: str, nested_paths: List[List[int]]):
        """ Delete some sectors of the nested path to conserve only the deepest (deeper mode), only
        the highest (higher mode) or both (bound mode)

        Parameters
        ----------
        mode: str
            Mode of path cutting
            - uncut: doesn't cut and keep all sectors
            - deeper: cut nested path and only conserve the deepest sector in the tree
            - higher: cut nested path and only conserve the highest sector in the tree
            - bound: cut nested path and only conserve the highest AND the deepest sectors in the
            tree
        nested_paths: List[List[int]]
            List of lists of nested path sectors indexes
        """
        to_del = []
        if mode == PATH_DEEPER:
            for path in nested_paths:
                to_del += path[:-1]
                to_keep = path[-1]
                root_p = self.parents[path[0]]
                self.parents[to_keep] = root_p
                self.labels[to_keep] = '... ' + self.labels[to_keep]
        elif mode == PATH_HIGHER:
            for path in nested_paths:
                to_del += path[1:]
                to_keep = path[0]
                to_keep_c = [i for i in range(self.len) if self.parents[i] == self.ids[path[-1]]]
                for c_i in to_keep_c:
                    self.parents[c_i] = self.ids[to_keep]
                self.labels[to_keep] += ' ...'
        elif mode == PATH_BOUND:
            for path in nested_paths:
                to_del += path[1:-1]
                to_keep_up = path[0]
                to_keep_do = path[-1]
                self.parents[to_keep_do] = self.ids[to_keep_up]
                if len(path) > 2:
                    self.labels[to_keep_up] += ' ...'
                    self.labels[to_keep_do] = '... ' + self.labels[to_keep_do]
        self.delete_value(to_del)

    def delete_value(self, v_index: int or List[int]):
        """ Delete a sector of TreeData from its index or a list of sectors from a list of indexes

        Parameters
        ----------
        v_index: int or List[int]
            Index or list of indexes of the sectors to delete
        """
        data = self.get_data_dict()
        if type(v_index) == int:
            v_index = [v_index]
        for i in sorted(v_index, reverse=True):
            for k, v in data.items():
                del v[i]
            self.len -= 1

    def get_col(self, index: int or List[int] = None) -> List or List[List]:
        """ Get a TreeData column from its index or a list of columns from a list of indexes.
        Column = all values associated with a sector.

        Parameters
        ----------
        index: int or List[int]
            Index or list of indexes of the sectors to get the column

        Returns
        -------
        List or List[List]
            Column or list of columns obtained from indexes
        """
        if index is None:
            index = list(range(self.len))
        if type(index) == int:
            index = [index]
        cols = list()
        for i in index:
            cols.append((self.ids[i], self.onto_ids[i], self.labels[i], self.parents[i],
                         self.count[i], self.ref_count[i], self.prop[i], self.ref_prop[i],
                         self.relative_prop[i], self.p_val[i]))
        return cols


# ==================================================================================================
# FUNCTIONS
# ==================================================================================================

def get_set2_abundance(set2_abundances: Dict[str, float] or None, c_label: str) -> float:
    """ Get the set2 abundance of a set1 concept.

    Parameters
    ----------
    set2_abundances: Dict[str, float] or None
        Dictionary associating to all the set1 concepts, its abundance
    c_label: str
        Label of the set2 concept

    Returns
    -------
    Abundance of the set1 concept. If set2 concept not in the set1, returns numpy.nan.
    """
    try:
        c_set2_abundance = set2_abundances[c_label]
    except KeyError:
        c_set2_abundance = np.nan
    return c_set2_abundance


def get_name(c_onto_id, names):
    if names is not None:
        try:
            c_label = names[c_onto_id]
        except KeyError:
            c_label = c_onto_id
    else:
        c_label = c_onto_id
    return c_label


def get_children_dict(parent_dict: Dict[str, List[str]]) -> Dict[str, List[str]]:
    """ Create the children dictionary from the parents dictionary.
    Parameters
    ----------
    parent_dict: Dict[str, List[str]]
        Dictionary associating for each class, its parents classes
    Returns
    -------
    Dict[str, List[str]]
        Dictionary associating for each class, its children classes
    """
    children_dict = dict()
    for c, ps in parent_dict.items():
        for p in ps:
            if p not in children_dict.keys():
                children_dict[p] = list()
            if c not in children_dict.keys():
                children_dict[c] = list()
            children_dict[p].append(c)
    return children_dict
