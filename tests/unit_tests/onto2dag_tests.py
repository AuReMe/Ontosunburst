import unittest
from unittest.mock import patch
import io
import sys
from functools import wraps
from ontosunburst.onto2dag import *

"""
Tests manually good file creation.
No automatic tests integrated.
"""

# ==================================================================================================
# GLOBAL
# ==================================================================================================

# --------------------------------------------------------------------------------------------------

CPT_LST = ['a', 'b', 'c']
RCPT_LST = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
CPT_AB = [1, 2, 3]
RCPT_AB = [1, 2, 3, 4, 5, 6, 7, 8]
ROOT = 'root'
ONTO_DAG = {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf'], 'd': ['cde'], 'e': ['cde', 'eg'],
           'f': ['cf'], 'g': ['gh', 'eg'], 'h': ['gh'],
           'ab': [ROOT], 'cde': ['cdecf', 'cdeeg'], 'cf': ['cdecf'],
           'eg': [ROOT, 'cdeeg'], 'gh': [ROOT],
           'cdecf': [ROOT], 'cdeeg': ['cdeeg+'], 'cdeeg+': [ROOT]}
ID2LAB = {ROOT: 'Root', 'cdeeg+': 'CDEEG+', 'cdeeg': 'CDEEG', 'cdecf': 'CDECF', 'gh': 'GH',
          'eg': 'EG', 'cde': 'CDE', 'cf': 'CF', 'h': 'H', 'g': 'G', 'f': 'F', 'e': 'E', 'd': 'D',
          'c': 'C', 'ab': 'AB', 'b': 'B'}

# ==================================================================================================
# FUNCTIONS UTILS
# ==================================================================================================
def dicts_with_sorted_lists_equal(dict1, dict2):
    if dict1.keys() != dict2.keys():
        return False
    for key in dict1:
        if sorted(dict1[key]) != sorted(dict2[key]):
            return False
    return True


def test_for(func):
    def decorator(test_func):
        @wraps(test_func)
        def wrapper(*args, **kwargs):
            return test_func(*args, **kwargs)

        wrapper._test_for = func
        return wrapper

    return decorator


class DualWriter(io.StringIO):
    def __init__(self, original_stdout):
        super().__init__()
        self.original_stdout = original_stdout

    def write(self, s):
        super().write(s)
        self.original_stdout.write(s)


# ==================================================================================================
# UNIT TESTS
# ==================================================================================================

# TESTS REDUCE DAG FUNCTIONS
# --------------------------------------------------------------------------------------------------
class TestReduceDAG(unittest.TestCase):

    @test_for(classify_concepts)
    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_classify_concepts_ok(self, mock_stdout):
        classified_concepts = classify_concepts(concepts=CPT_LST, ontology_dag=ONTO_DAG)
        output = mock_stdout.getvalue().strip()
        self.assertEqual(output, '3 concepts to classify\n'
                                 '3/3 concepts classified')
        self.assertEqual(classified_concepts, {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf']})



    @test_for(classify_concepts)
    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_classify_concepts_errors(self, mock_stdout):
        classified_concepts = classify_concepts(concepts=CPT_LST + ['x'], ontology_dag=ONTO_DAG)
        output = mock_stdout.getvalue().strip()
        self.assertEqual(output, '4 concepts to classify\n'
                                 'x not classified.\n'
                                 '3/4 concepts classified')
        self.assertEqual(classified_concepts, {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf']})

    @test_for(get_parents)
    def test_get_parents_linear_path(self):
        # Simple linear direction
        parents = get_parents('a', {'ab'}, ONTO_DAG, ROOT)
        self.assertEqual(parents, {'root', 'ab'})

    @test_for(get_parents)
    def test_get_parents_complex_path(self):
        # With multiple parents having multiple parents and different size of path until root
        parents = get_parents('c', {'cde', 'cf'}, ONTO_DAG, ROOT)
        self.assertEqual(parents, {'cdeeg+', 'root', 'cf', 'cde', 'cdecf', 'cdeeg'})

    @test_for(get_all_classes)
    def test_get_all_classes(self):
        leaf_classes = {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf']}
        all_classes_met = get_all_classes(leaf_classes, ONTO_DAG, ROOT)
        wanted_all_classes = {'a': {'root', 'ab'}, 'b': {'root', 'ab'},
                              'c': {'cdeeg+', 'cde', 'cdeeg', 'root', 'cdecf', 'cf'}}
        self.assertEqual(all_classes_met, wanted_all_classes)


# TESTS WEIGHTS CALCULATION
# --------------------------------------------------------------------------------------------------
class TestWeightsCalculation(unittest.TestCase):
    @test_for(get_abundance_dict)
    def test_get_abundance_dict_abundances(self):
        abundance_dict = get_abundance_dict(abundances=CPT_AB, concepts=CPT_LST)
        self.assertEqual(abundance_dict, {'a': 1, 'b': 2, 'c': 3})

    @test_for(get_abundance_dict)
    def test_get_abundance_dict_no_abundances(self):
        abundance_dict = get_abundance_dict(abundances=None, concepts=CPT_LST)
        self.assertEqual(abundance_dict, {'a': 1, 'b': 1, 'c': 1})

    @test_for(get_abundance_dict)
    def test_get_abundance_dict_abundances_ref(self):
        abundance_dict = get_abundance_dict(abundances=RCPT_AB, concepts=RCPT_LST)
        self.assertEqual(abundance_dict, {'a': 1, 'b': 2, 'c': 3, 'd': 4,
                                          'e': 5, 'f': 6, 'g': 7, 'h': 8})

    @test_for(get_abundance_dict)
    def test_get_abundance_dict_no_abundances_ref(self):
        abundance_dict = get_abundance_dict(abundances=None, concepts=RCPT_LST)
        self.assertEqual(abundance_dict, {'a': 1, 'b': 1, 'c': 1, 'd': 1,
                                          'e': 1, 'f': 1, 'g': 1, 'h': 1})

    @test_for(get_abundance_dict)
    def test_get_abundance_dict_errors(self):
        with self.assertRaises(AttributeError) as e:
            get_abundance_dict(abundances=CPT_AB + [4], concepts=CPT_LST)
        self.assertEqual(str(e.exception), 'Length of concepts IDs list must '
                                           'be equal to its abundances list length : 3 != 4')


    @test_for(calculate_weights)
    def test_get_classes_abundance_leaves(self):
        all_classes = {'a': {'root', 'ab'}, 'b': {'root', 'ab'},
                       'c': {'cdecf', 'cdeeg+', 'root', 'cde', 'cdeeg', 'cf'},
                       'd': {'cdecf', 'cdeeg+', 'root', 'cde', 'cdeeg'},
                       'e': {'cdeeg+', 'root', 'cde', 'cdecf', 'eg', 'cdeeg'},
                       'f': {'cdecf', 'root', 'cf'},
                       'g': {'cdeeg', 'cdeeg+', 'root', 'eg', 'gh'}, 'h': {'root', 'gh'}}
        abundances_dict = {'a': 1, 'b': 2, 'c': 3, 'd': 4, 'e': 5, 'f': 6, 'g': 7, 'h': 8}
        classes_abundances = calculate_weights(all_classes, abundances_dict, show_leaves=True)
        wanted_abundances = {'root': 36, 'cdeeg+': 19, 'cdeeg': 19, 'cdecf': 18, 'gh': 15,
                             'eg': 12, 'cde': 12, 'cf': 9, 'h': 8, 'g': 7, 'f': 6, 'e': 5,
                             'd': 4, 'c': 3, 'ab': 3, 'b': 2, 'a': 1}
        self.assertEqual(classes_abundances, wanted_abundances)

    @test_for(calculate_weights)
    def test_get_classes_abundance_leaves_sub(self):
        all_classes = {'a': {'root', 'ab'}, 'b': {'root', 'ab'},
                       'c': {'cdeeg+', 'cde', 'cdeeg', 'root', 'cdecf', 'cf'}}
        abundances_dict = {'a': 1, 'b': 2, 'c': 3}
        classes_abundances = calculate_weights(all_classes, abundances_dict, show_leaves=True)
        wanted_abundances = {'root': 6, 'cde': 3, 'cf': 3, 'cdecf': 3, 'cdeeg+': 3, 'cdeeg': 3,
                             'c': 3, 'ab': 3, 'b': 2, 'a': 1}
        self.assertEqual(classes_abundances, wanted_abundances)

    @test_for(calculate_weights)
    def test_get_classes_abundance_no_leaves(self):
        all_classes = {'a': {'root', 'ab'}, 'b': {'root', 'ab'},
                       'c': {'cdecf', 'cdeeg+', 'root', 'cde', 'cdeeg', 'cf'},
                       'd': {'cdecf', 'cdeeg+', 'root', 'cde', 'cdeeg'},
                       'e': {'cdeeg+', 'root', 'cde', 'cdecf', 'eg', 'cdeeg'},
                       'f': {'cdecf', 'root', 'cf'},
                       'g': {'cdeeg', 'cdeeg+', 'root', 'eg', 'gh'}, 'h': {'root', 'gh'}}
        abundances_dict = {'a': 1, 'b': 2, 'c': 3, 'd': 4, 'e': 5, 'f': 6, 'g': 7, 'h': 8}
        classes_abundances = calculate_weights(all_classes, abundances_dict, show_leaves=False)
        wanted_abundances = {'root': 36, 'cdeeg+': 19, 'cdeeg': 19, 'cdecf': 18, 'gh': 15,
                             'eg': 12, 'cde': 12, 'cf': 9, 'ab': 3}
        self.assertEqual(classes_abundances, wanted_abundances)

    @test_for(calculate_weights)
    def test_get_classes_abundance_no_leaves_sub(self):
        all_classes = {'a': {'root', 'ab'}, 'b': {'root', 'ab'},
                       'c': {'cdeeg+', 'cde', 'cdeeg', 'root', 'cdecf', 'cf'}}
        abundances_dict = {'a': 1, 'b': 2, 'c': 3}
        classes_abundances = calculate_weights(all_classes, abundances_dict, show_leaves=False)
        wanted_abundances = {'root': 6, 'cde': 3, 'cf': 3, 'cdecf': 3, 'cdeeg+': 3, 'cdeeg': 3,
                             'ab': 3}
        self.assertEqual(classes_abundances, wanted_abundances)

    @test_for(calculate_weights)
    def test_get_classes_abundance_different_level_abundances(self):
        all_classes = {'c': {'cdecf', 'cdeeg+', 'root', 'cde', 'cdeeg', 'cf'},
                       'd': {'cdecf', 'cdeeg+', 'root', 'cde', 'cdeeg'},
                       'e': {'cdeeg+', 'root', 'cde', 'cdecf', 'eg', 'cdeeg'},
                       'f': {'cdecf', 'root', 'cf'}, 'cf': {'cdecf', 'root'}}
        abundances_dict = {'c': 3, 'd': 4, 'e': 5, 'f': 2, 'cf': 2}
        classes_abundances = calculate_weights(all_classes, abundances_dict, show_leaves=True)
        wanted_abundances = {ROOT: 16, 'cdecf': 16, 'cde': 12, 'cdeeg+': 12, 'cdeeg': 12,
                             'cf': 7, 'eg': 5, 'e': 5, 'd': 4, 'c': 3, 'f': 2}
        self.assertEqual(classes_abundances, wanted_abundances)


# TEST MAIN FUNCTIONS
# --------------------------------------------------------------------------------------------------

class TestMainFunctions(unittest.TestCase):

    @test_for(reduce_d_ontology)
    def test_reduce_d_ontology_onto_dag(self):
        classes_abundance = {'root': 6, 'cde': 3, 'cf': 3, 'cdecf': 3, 'cdeeg+': 3, 'cdeeg': 3,
                             'ab': 3}
        d_ontology_reduced = reduce_d_ontology(ONTO_DAG, classes_abundance)
        wanted_d_ontology_reduced = {'ab': ['root'], 'cde': ['cdecf', 'cdeeg'],
                                     'cf': ['cdecf'], 'cdecf': ['root'], 'cdeeg': ['cdeeg+'],
                                     'cdeeg+': ['root']}
        self.assertEqual(d_ontology_reduced, wanted_d_ontology_reduced)

    @test_for(reduce_d_ontology)
    def test_reduce_d_ontology_id_to_label(self):
        classes_abundance = {'root': 6, 'cde': 3, 'cf': 3, 'cdecf': 3, 'cdeeg+': 3, 'cdeeg': 3,
                             'ab': 3}
        d_ontology_reduced = reduce_d_ontology(ID2LAB, classes_abundance)
        wanted_d_ontology_reduced = {'root': 'Root', 'cdeeg+': 'CDEEG+', 'cdeeg': 'CDEEG',
                                     'cdecf': 'CDECF', 'cde': 'CDE', 'cf': 'CF', 'ab': 'AB'}
        self.assertEqual(d_ontology_reduced, wanted_d_ontology_reduced)

    @test_for(ontology_to_weighted_dag)
    def test_ontology_to_weighted_dag_no_lvs(self):
        calculated_weights = ontology_to_weighted_dag(CPT_LST, CPT_AB, ROOT, ONTO_DAG, False)
        wanted_abundances = {'root': 6, 'cde': 3, 'cf': 3, 'cdecf': 3, 'cdeeg+': 3, 'cdeeg': 3,
                             'ab': 3}
        self.assertEqual(calculated_weights, wanted_abundances)

    @test_for(ontology_to_weighted_dag)
    def test_ontology_to_weighted_dag_ref_lvs(self):
        calculated_weights = ontology_to_weighted_dag(RCPT_LST, RCPT_AB, ROOT, ONTO_DAG, True)
        wanted_abundances = {'root': 36, 'cdeeg+': 19, 'cdeeg': 19, 'cdecf': 18, 'gh': 15,
                             'eg': 12, 'cde': 12, 'cf': 9, 'h': 8, 'g': 7, 'f': 6, 'e': 5,
                             'd': 4, 'c': 3, 'ab': 3, 'b': 2, 'a': 1}
        self.assertEqual(calculated_weights, wanted_abundances)
