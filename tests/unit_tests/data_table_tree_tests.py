import copy
import unittest
import io

from numpy import nan
from functools import wraps
from ontosunburst.data_table_tree import *
from ontosunburst.ontology import *

"""
Tests manually good file creation.
No automatic tests integrated.
"""

# ==================================================================================================
# GLOBAL
# ==================================================================================================

# --------------------------------------------------------------------------------------------------

MC_ONTO = {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf'], 'd': ['cde'], 'e': ['cde', 'eg'],
           'f': ['cf'], 'g': ['gh', 'eg'], 'h': ['gh'],
           'ab': [ROOTS[METACYC]], 'cde': ['cdecf', 'cdeeg'], 'cf': ['cdecf'],
           'eg': ['cdeeg', ROOTS[METACYC]], 'gh': [ROOTS[METACYC]],
           'cdecf': [ROOTS[METACYC]], 'cdeeg': ['cdeeg+'], 'cdeeg+': [ROOTS[METACYC]]}

MET_LAB_D = {'FRAMES': 6, 'cde': 3, 'cf': 3, 'cdecf': 3, 'cdeeg+': 3, 'cdeeg': 3, 'c': 3, 'ab': 3,
             'b': 2, 'a': 1}
MET_RAB_D = {'FRAMES': 36, 'cdeeg+': 19, 'cdeeg': 19, 'cdecf': 18, 'gh': 15, 'eg': 12, 'cde': 12,
             'cf': 9, 'h': 8, 'g': 7, 'f': 6, 'e': 5, 'd': 4, 'c': 3, 'ab': 3, 'b': 2, 'a': 1}

NAMES = {'FRAMES': 'Root', 'cdeeg+': 'CDEEG+', 'cdeeg': 'CDEEG', 'cdecf': 'CDECF', 'gh': 'GH',
         'eg': 'EG', 'cde': 'CDE', 'cf': 'CF', 'h': 'H', 'g': 'G', 'f': 'F', 'e': 'E', 'd': 'D',
         'c': 'C', 'ab': 'AB', 'b': 'B'}

DATA = {'ID': ['FRAMES', 'cdeeg+__FRAMES', 'cdeeg__cdeeg+__FRAMES', 'cdecf__FRAMES', 'gh__FRAMES',
               'eg__cdeeg__cdeeg+__FRAMES', 'eg__FRAMES', 'cde__cdeeg__cdeeg+__FRAMES',
               'cde__cdecf__FRAMES', 'cf__cdecf__FRAMES', 'h__gh__FRAMES',
               'g__eg__cdeeg__cdeeg+__FRAMES', 'g__eg__FRAMES', 'g__gh__FRAMES',
               'f__cf__cdecf__FRAMES', 'e__eg__cdeeg__cdeeg+__FRAMES',
               'e__cde__cdeeg__cdeeg+__FRAMES', 'e__eg__FRAMES', 'e__cde__cdecf__FRAMES',
               'd__cde__cdecf__FRAMES', 'd__cde__cdeeg__cdeeg+__FRAMES', 'c__cde__cdecf__FRAMES',
               'c__cf__cdecf__FRAMES', 'c__cde__cdeeg__cdeeg+__FRAMES', 'ab__FRAMES',
               'b__ab__FRAMES', 'a__ab__FRAMES'],
        'Onto ID': ['FRAMES', 'cdeeg+', 'cdeeg', 'cdecf', 'gh', 'eg', 'eg', 'cde', 'cde', 'cf', 'h',
                    'g', 'g', 'g', 'f', 'e', 'e', 'e', 'e', 'd', 'd', 'c', 'c', 'c', 'ab', 'b',
                    'a'],
        'Parent': ['', 'FRAMES', 'cdeeg+__FRAMES', 'FRAMES', 'FRAMES', 'cdeeg__cdeeg+__FRAMES',
                   'FRAMES', 'cdeeg__cdeeg+__FRAMES', 'cdecf__FRAMES', 'cdecf__FRAMES',
                   'gh__FRAMES', 'eg__cdeeg__cdeeg+__FRAMES', 'eg__FRAMES', 'gh__FRAMES',
                   'cf__cdecf__FRAMES', 'eg__cdeeg__cdeeg+__FRAMES', 'cde__cdeeg__cdeeg+__FRAMES',
                   'eg__FRAMES', 'cde__cdecf__FRAMES', 'cde__cdecf__FRAMES',
                   'cde__cdeeg__cdeeg+__FRAMES', 'cde__cdecf__FRAMES', 'cf__cdecf__FRAMES',
                   'cde__cdeeg__cdeeg+__FRAMES', 'FRAMES', 'ab__FRAMES', 'ab__FRAMES'],
        'Label': ['FRAMES', 'cdeeg+', 'cdeeg', 'cdecf', 'gh', 'eg', 'eg', 'cde', 'cde', 'cf', 'h',
                  'g', 'g', 'g', 'f', 'e', 'e', 'e', 'e', 'd', 'd', 'c', 'c', 'c', 'ab', 'b', 'a'],
        'Count': [6, 3, 3, 3, nan, nan, nan, 3, 3, 3, nan, nan, nan, nan, nan, nan, nan, nan, nan,
                  nan, nan, 3, 3, 3, 3, 2, 1],
        'Reference count': [36, 19, 19, 18, 15, 12, 12, 12, 12, 9, 8, 7, 7, 7, 6, 5, 5, 5, 5, 4, 4,
                            3, 3, 3, 3, 2, 1]}

W_PROP = [1.0, 0.5, 0.5, 0.5, nan, nan, nan, 0.5, 0.5, 0.5, nan, nan, nan, nan, nan, nan, nan, nan,
          nan, nan, nan, 0.5, 0.5, 0.5, 0.5, 0.3333333333333333, 0.16666666666666666]

W_REF_PROP = [1.0, 0.5277777777777778, 0.5277777777777778, 0.5, 0.4166666666666667,
              0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.25,
              0.2222222222222222, 0.19444444444444445, 0.19444444444444445, 0.19444444444444445,
              0.16666666666666666, 0.1388888888888889, 0.1388888888888889, 0.1388888888888889,
              0.1388888888888889, 0.1111111111111111, 0.1111111111111111, 0.08333333333333333,
              0.08333333333333333, 0.08333333333333333, 0.08333333333333333, 0.05555555555555555,
              0.027777777777777776]

W_RELAT_PROP = [1000000, 283582, 283582, 268656, 223880, 141791, 179104, 141791, 153517, 115138,
                119402, 82711, 104477, 104477, 76758, 59079, 59079, 74626, 63965, 51172, 47263,
                38379, 38379, 35447, 44776, 29850, 14925]
W_REL_PROP = {'FRAMES': 1000000, 'cdeeg+__FRAMES': 283582, 'cdeeg__cdeeg+__FRAMES': 283582,
              'cdecf__FRAMES': 268656, 'gh__FRAMES': 223880, 'eg__cdeeg__cdeeg+__FRAMES': 141791,
              'eg__FRAMES': 179104, 'cde__cdeeg__cdeeg+__FRAMES': 141791,
              'cde__cdecf__FRAMES': 153517, 'cf__cdecf__FRAMES': 115138, 'h__gh__FRAMES': 119402,
              'g__eg__cdeeg__cdeeg+__FRAMES': 82711, 'g__eg__FRAMES': 104477,
              'g__gh__FRAMES': 104477, 'f__cf__cdecf__FRAMES': 76758,
              'e__eg__cdeeg__cdeeg+__FRAMES': 59079, 'e__cde__cdeeg__cdeeg+__FRAMES': 59079,
              'e__eg__FRAMES': 74626, 'e__cde__cdecf__FRAMES': 63965,
              'd__cde__cdecf__FRAMES': 51172, 'd__cde__cdeeg__cdeeg+__FRAMES': 47263,
              'c__cde__cdecf__FRAMES': 38379, 'c__cf__cdecf__FRAMES': 38379,
              'c__cde__cdeeg__cdeeg+__FRAMES': 35447, 'ab__FRAMES': 44776, 'b__ab__FRAMES': 29850,
              'a__ab__FRAMES': 14925}

DATA_PROP = DATA
DATA_PROP[PROP] = W_PROP
DATA_PROP[REF_PROP] = W_REF_PROP

DATA_R_PROP = DATA_PROP
DATA_R_PROP[RELAT_PROP] = W_RELAT_PROP


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


def data_to_lines(dico):
    lines = set()
    for i in range(len(dico[IDS])):
        line = (dico[IDS][i], dico[ONTO_ID][i], dico[PARENT][i], dico[LABEL][i], dico[COUNT][i],
                dico[REF_COUNT][i])
        if PROP in dico:
            line = line + (dico[PROP][i],)
        if REF_PROP in dico:
            line = line + (dico[REF_PROP][i],)
        if RELAT_PROP in dico:
            line = line + (dico[RELAT_PROP][i],)
        if PVAL in dico:
            line = line + (dico[PVAL][i],)
        lines.add(line)
    return lines


# ==================================================================================================
# UNIT TESTS
# ==================================================================================================

# TEST
# --------------------------------------------------------------------------------------------------

class TestGenerateDataTable(unittest.TestCase):

    @test_for(get_sub_abundance)
    def test_get_sub_abundances_exists_diff(self):
        sub_abu = get_sub_abundance(MET_LAB_D, 'cf', 9)
        self.assertEqual(sub_abu, 3)

    @test_for(get_sub_abundance)
    def test_get_sub_abundances_exists_equ(self):
        sub_abu = get_sub_abundance(MET_LAB_D, 'a', 1)
        self.assertEqual(sub_abu, 1)

    @test_for(get_sub_abundance)
    def test_get_sub_abundances_not_exists(self):
        sub_abu = get_sub_abundance(MET_LAB_D, 'eg', 12)
        self.assertTrue(np.isnan(sub_abu))

    @test_for(get_sub_abundance)
    def test_get_sub_abundances_no_sub(self):
        sub_abu = get_sub_abundance(None, 'eg', 12)
        self.assertEqual(sub_abu, 12)

    @test_for(add_value_data)
    def test_add_value_data(self):
        data = DataTable()
        data.add_value(m_id='bjr', onto_id='Bjr_0', label='bonjour', count=2, ref_count=8,
                       parent='salutations')
        data.add_value(m_id='slt', onto_id='sl_1', label='salut', count=0.5, ref_count=2.3,
                       parent='salutations')
        wanted_data = {IDS: ['bjr', 'slt'],
                       ONTO_ID: ['Bjr_0', 'sl_1'],
                       PARENT: ['salutations', 'salutations'],
                       LABEL: ['bonjour', 'salut'],
                       COUNT: [2, 0.5],
                       REF_COUNT: [8, 2.3],
                       PROP: [nan, nan], REF_PROP: [nan, nan], RELAT_PROP: [nan, nan],
                       PVAL: [nan, nan]}
        self.assertEqual(data.get_data_dict(), wanted_data)

    @test_for(get_all_ids)
    def test_get_all_c_ids(self):
        all_ids = get_all_ids('c', 'c', MC_ONTO, ROOTS[METACYC], set())
        wanted_ids = {'c__cf__cdecf__FRAMES', 'c__cde__cdecf__FRAMES',
                      'c__cde__cdeeg__cdeeg+__FRAMES'}
        self.assertEqual(all_ids, wanted_ids)

    @test_for(get_all_ids)
    def test_get_all_e_ids(self):
        all_ids = get_all_ids('e', 'e', MC_ONTO, ROOTS[METACYC], set())
        wanted_ids = {'e__cde__cdeeg__cdeeg+__FRAMES', 'e__eg__FRAMES',
                      'e__eg__cdeeg__cdeeg+__FRAMES', 'e__cde__cdecf__FRAMES'}
        self.assertEqual(all_ids, wanted_ids)

    @test_for(get_all_ids)
    def test_get_all_eg_ids(self):
        all_ids = get_all_ids('eg', 'eg', MC_ONTO, ROOTS[METACYC], set())
        wanted_ids = {'eg__FRAMES', 'eg__cdeeg__cdeeg+__FRAMES'}
        self.assertEqual(all_ids, wanted_ids)

    @test_for(get_fig_parameters)
    def test_get_fig_parameters(self):
        data = DataTable()
        data.fill_parameters(classes_abondance=MET_RAB_D, parent_dict=MC_ONTO,
                             root_item=ROOTS[METACYC], subset_abundance=MET_LAB_D, names=None)
        lines = set(data.get_col())
        w_lines = {('cde__cdeeg__cdeeg+__FRAMES', 'cde', 'cde', 'cdeeg__cdeeg+__FRAMES', 3, 12, nan, nan, nan, nan), ('d__cde__cdeeg__cdeeg+__FRAMES', 'd', 'd', 'cde__cdeeg__cdeeg+__FRAMES', nan, 4, nan, nan, nan, nan), ('b__ab__FRAMES', 'b', 'b', 'ab__FRAMES', 2, 2, nan, nan, nan, nan), ('g__eg__FRAMES', 'g', 'g', 'eg__FRAMES', nan, 7, nan, nan, nan, nan), ('eg__FRAMES', 'eg', 'eg', 'FRAMES', nan, 12, nan, nan, nan, nan), ('c__cde__cdecf__FRAMES', 'c', 'c', 'cde__cdecf__FRAMES', 3, 3, nan, nan, nan, nan), ('g__eg__cdeeg__cdeeg+__FRAMES', 'g', 'g', 'eg__cdeeg__cdeeg+__FRAMES', nan, 7, nan, nan, nan, nan), ('e__eg__cdeeg__cdeeg+__FRAMES', 'e', 'e', 'eg__cdeeg__cdeeg+__FRAMES', nan, 5, nan, nan, nan, nan), ('cdeeg__cdeeg+__FRAMES', 'cdeeg', 'cdeeg', 'cdeeg+__FRAMES', 3, 19, nan, nan, nan, nan), ('cdecf__FRAMES', 'cdecf', 'cdecf', 'FRAMES', 3, 18, nan, nan, nan, nan), ('ab__FRAMES', 'ab', 'ab', 'FRAMES', 3, 3, nan, nan, nan, nan), ('e__cde__cdeeg__cdeeg+__FRAMES', 'e', 'e', 'cde__cdeeg__cdeeg+__FRAMES', nan, 5, nan, nan, nan, nan), ('c__cde__cdeeg__cdeeg+__FRAMES', 'c', 'c', 'cde__cdeeg__cdeeg+__FRAMES', 3, 3, nan, nan, nan, nan), ('c__cf__cdecf__FRAMES', 'c', 'c', 'cf__cdecf__FRAMES', 3, 3, nan, nan, nan, nan), ('cde__cdecf__FRAMES', 'cde', 'cde', 'cdecf__FRAMES', 3, 12, nan, nan, nan, nan), ('FRAMES', 'FRAMES', 'FRAMES', '', 6, 36, nan, nan, nan, nan), ('f__cf__cdecf__FRAMES', 'f', 'f', 'cf__cdecf__FRAMES', nan, 6, nan, nan, nan, nan), ('eg__cdeeg__cdeeg+__FRAMES', 'eg', 'eg', 'cdeeg__cdeeg+__FRAMES', nan, 12, nan, nan, nan, nan), ('e__eg__FRAMES', 'e', 'e', 'eg__FRAMES', nan, 5, nan, nan, nan, nan), ('g__gh__FRAMES', 'g', 'g', 'gh__FRAMES', nan, 7, nan, nan, nan, nan), ('h__gh__FRAMES', 'h', 'h', 'gh__FRAMES', nan, 8, nan, nan, nan, nan), ('cdeeg+__FRAMES', 'cdeeg+', 'cdeeg+', 'FRAMES', 3, 19, nan, nan, nan, nan), ('gh__FRAMES', 'gh', 'gh', 'FRAMES', nan, 15, nan, nan, nan, nan), ('e__cde__cdecf__FRAMES', 'e', 'e', 'cde__cdecf__FRAMES', nan, 5, nan, nan, nan, nan), ('cf__cdecf__FRAMES', 'cf', 'cf', 'cdecf__FRAMES', 3, 9, nan, nan, nan, nan), ('a__ab__FRAMES', 'a', 'a', 'ab__FRAMES', 1, 1, nan, nan, nan, nan), ('d__cde__cdecf__FRAMES', 'd', 'd', 'cde__cdecf__FRAMES', nan, 4, nan, nan, nan, nan)}
        self.assertEqual(lines, w_lines)

    @test_for(get_fig_parameters)
    def test_get_fig_parameters_names(self):
        data = DataTable()
        data.fill_parameters(classes_abondance=MET_RAB_D, parent_dict=MC_ONTO,
                             root_item=ROOTS[METACYC], subset_abundance=MET_LAB_D, names=NAMES)
        lines = set(data.get_col())
        w_lines = {('ab__FRAMES', 'ab', 'AB', 'FRAMES', 3, 3, nan, nan, nan, nan), ('d__cde__cdecf__FRAMES', 'd', 'D', 'cde__cdecf__FRAMES', nan, 4, nan, nan, nan, nan), ('b__ab__FRAMES', 'b', 'B', 'ab__FRAMES', 2, 2, nan, nan, nan, nan), ('h__gh__FRAMES', 'h', 'H', 'gh__FRAMES', nan, 8, nan, nan, nan, nan), ('g__eg__cdeeg__cdeeg+__FRAMES', 'g', 'G', 'eg__cdeeg__cdeeg+__FRAMES', nan, 7, nan, nan, nan, nan), ('cdeeg+__FRAMES', 'cdeeg+', 'CDEEG+', 'FRAMES', 3, 19, nan, nan, nan, nan), ('FRAMES', 'FRAMES', 'FRAMES', '', 6, 36, nan, nan, nan, nan), ('g__eg__FRAMES', 'g', 'G', 'eg__FRAMES', nan, 7, nan, nan, nan, nan), ('cde__cdeeg__cdeeg+__FRAMES', 'cde', 'CDE', 'cdeeg__cdeeg+__FRAMES', 3, 12, nan, nan, nan, nan), ('eg__cdeeg__cdeeg+__FRAMES', 'eg', 'EG', 'cdeeg__cdeeg+__FRAMES', nan, 12, nan, nan, nan, nan), ('c__cde__cdeeg__cdeeg+__FRAMES', 'c', 'C', 'cde__cdeeg__cdeeg+__FRAMES', 3, 3, nan, nan, nan, nan), ('g__gh__FRAMES', 'g', 'G', 'gh__FRAMES', nan, 7, nan, nan, nan, nan), ('e__eg__FRAMES', 'e', 'E', 'eg__FRAMES', nan, 5, nan, nan, nan, nan), ('eg__FRAMES', 'eg', 'EG', 'FRAMES', nan, 12, nan, nan, nan, nan), ('e__eg__cdeeg__cdeeg+__FRAMES', 'e', 'E', 'eg__cdeeg__cdeeg+__FRAMES', nan, 5, nan, nan, nan, nan), ('c__cf__cdecf__FRAMES', 'c', 'C', 'cf__cdecf__FRAMES', 3, 3, nan, nan, nan, nan), ('gh__FRAMES', 'gh', 'GH', 'FRAMES', nan, 15, nan, nan, nan, nan), ('cdeeg__cdeeg+__FRAMES', 'cdeeg', 'CDEEG', 'cdeeg+__FRAMES', 3, 19, nan, nan, nan, nan), ('e__cde__cdeeg__cdeeg+__FRAMES', 'e', 'E', 'cde__cdeeg__cdeeg+__FRAMES', nan, 5, nan, nan, nan, nan), ('d__cde__cdeeg__cdeeg+__FRAMES', 'd', 'D', 'cde__cdeeg__cdeeg+__FRAMES', nan, 4, nan, nan, nan, nan), ('c__cde__cdecf__FRAMES', 'c', 'C', 'cde__cdecf__FRAMES', 3, 3, nan, nan, nan, nan), ('a__ab__FRAMES', 'a', 'a', 'ab__FRAMES', 1, 1, nan, nan, nan, nan), ('e__cde__cdecf__FRAMES', 'e', 'E', 'cde__cdecf__FRAMES', nan, 5, nan, nan, nan, nan), ('cdecf__FRAMES', 'cdecf', 'CDECF', 'FRAMES', 3, 18, nan, nan, nan, nan), ('cde__cdecf__FRAMES', 'cde', 'CDE', 'cdecf__FRAMES', 3, 12, nan, nan, nan, nan), ('cf__cdecf__FRAMES', 'cf', 'CF', 'cdecf__FRAMES', 3, 9, nan, nan, nan, nan), ('f__cf__cdecf__FRAMES', 'f', 'F', 'cf__cdecf__FRAMES', nan, 6, nan, nan, nan, nan)}
        self.assertEqual(lines, w_lines)


class TestAddProportionDataTable(unittest.TestCase):

    @test_for(get_data_proportion)
    def test_get_data_proportion_no_relative(self):
        data = DataTable()
        data.fill_parameters(classes_abondance=MET_RAB_D, parent_dict=MC_ONTO,
                             root_item=ROOTS[METACYC], subset_abundance=MET_LAB_D, names=NAMES)
        data.calculate_proportions(False)
        for i in range(data.len):
            if np.isnan(data.prop[i]):
                self.assertTrue(np.isnan(W_PROP[i]))
            else:
                self.assertEqual(data.prop[i], W_PROP[i])

    @test_for(get_data_proportion)
    def test_get_data_proportion_no_relative_ref(self):
        data = DataTable()
        data.fill_parameters(classes_abondance=MET_RAB_D, parent_dict=MC_ONTO,
                             root_item=ROOTS[METACYC], subset_abundance=MET_LAB_D, names=NAMES)
        data.calculate_proportions(False)
        for i in range(data.len):
            self.assertEqual(data.ref_prop[i], W_REF_PROP[i])

    @test_for(get_data_proportion)
    def test_get_data_proportion_relative(self):
        data = DataTable()
        data.fill_parameters(classes_abondance=MET_RAB_D, parent_dict=MC_ONTO,
                             root_item=ROOTS[METACYC], subset_abundance=MET_LAB_D, names=NAMES)
        data.calculate_proportions(True)
        for k, v in W_REL_PROP.items():
            self.assertEqual(data.relative_prop[data.ids.index(k)], v)


# ENRICHMENT TESTS
# ==================================================================================================

ENRICH_DATA = {IDS: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
               ONTO_ID: ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09'],
               PARENT: ['', 0, 0, 0, 0, 1, 1, 1, 2, 2],
               LABEL: ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09'],
               COUNT: [50, 5, 25, 20, 1, 5, nan, nan, 1, 1],
               REF_COUNT: [100, 40, 30, 20, 10, 20, 5, 1, 1, 3],
               PROP: [1, 0.1, 0.5, 0.4, 0.02, 0.1, nan, nan, 0.02, 0.02],
               REF_PROP: [1, 0.4, 0.3, 0.2, 0.1, 0.2, 0.05, 0.01, 0.01, 0.03],
               RELAT_PROP: [1, 0.4, 0.3, 0.2, 0.1, 0.2, 0.05, 0.01, 0.01, 0.03]}
ENRICH_REF_AB = {'00': 100, '01': 40, '02': 30, '03': 20, '04': 10, '05': 20, '06': 5, '07': 1,
                 '08': 1, '09': 3}
LABEL_NAMES = ['r', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine']


# Expected :
# Over : 2, 3 | Under : 1, 4, 5 | No diff : 0, 8, 9 | Nan : 6, 7


class TestEnrichmentAnalysis(unittest.TestCase):

    @test_for(get_data_enrichment_analysis)
    def test_get_data_enrichment_analysis_single_value(self):
        data, significant = get_data_enrichment_analysis(ENRICH_DATA, ENRICH_REF_AB, BINOMIAL_TEST)
        lines = data_to_lines(data)
        p_value_1 = [l[-1] for l in lines if l[0] == 1][0]
        M = 100
        N = 50
        m = 40
        n = 5
        exp_p_value_1 = stats.binomtest(n, N, m / M, alternative='two-sided').pvalue
        exp_p_value_1 = np.log10(exp_p_value_1)
        self.assertEqual(p_value_1, exp_p_value_1)

    @test_for(get_data_enrichment_analysis)
    def test_get_data_enrichment_analysis_binomial(self):
        data, significant = get_data_enrichment_analysis(ENRICH_DATA, ENRICH_REF_AB, BINOMIAL_TEST)
        lines = data_to_lines(data)
        exp_significant = {'01': 3.7996e-06, '03': 0.0011251149, '02': 0.0030924096}
        exp_lines = {(7, '07', 1, '07', nan, 1, nan, 0.01, 0.01, nan),
                     (3, '03', 0, '03', 20, 20, 0.4, 0.2, 0.2, 2.948803113091024),
                     (5, '05', 1, '05', 5, 20, 0.1, 0.2, 0.2, -1.103304935668835),
                     (4, '04', 0, '04', 1, 10, 0.02, 0.1, 0.1, -1.2341542222355069),
                     (2, '02', 0, '02', 25, 30, 0.5, 0.3, 0.3, 2.509702991379166),
                     (1, '01', 0, '01', 5, 40, 0.1, 0.4, 0.4, -5.420266413988895),
                     (9, '09', 2, '09', 1, 3, 0.02, 0.03, 0.03, 0.0),
                     (8, '08', 2, '08', 1, 1, 0.02, 0.01, 0.01, 0.4034095751193356),
                     (0, '00', '', '00', 50, 100, 1, 1, 1, 0.0),
                     (6, '06', 1, '06', nan, 5, nan, 0.05, 0.05, nan)}
        self.assertEqual(lines, exp_lines)
        self.assertEqual(significant, exp_significant)

    @test_for(get_data_enrichment_analysis)
    def test_get_data_enrichment_analysis_hypergeometric(self):
        data, significant = get_data_enrichment_analysis(ENRICH_DATA, ENRICH_REF_AB, HYPERGEO_TEST)
        lines = data_to_lines(data)
        exp_lines = {(8, '08', 2, '08', 1, 1, 0.02, 0.01, 0.01, -0.0),
                     (3, '03', 0, '03', 20, 20, 0.4, 0.2, 0.2, 6.754831139005899),
                     (6, '06', 1, '06', nan, 5, nan, 0.05, 0.05, nan),
                     (5, '05', 1, '05', 5, 20, 0.1, 0.2, 0.2, -1.6413993451973743),
                     (9, '09', 2, '09', 1, 3, 0.02, 0.03, 0.03, -1.4464911998299308e-16),
                     (4, '04', 0, '04', 1, 10, 0.02, 0.1, 0.1, -1.8051946563380086),
                     (2, '02', 0, '02', 25, 30, 0.5, 0.3, 0.3, 4.692610428021241),
                     (0, '00', '', '00', 50, 100, 1, 1, 1, 0.3010299956639812),
                     (1, '01', 0, '01', 5, 40, 0.1, 0.4, 0.4, -9.138873998573988),
                     (7, '07', 1, '07', nan, 1, nan, 0.01, 0.01, nan)}
        exp_significant = {'01': 7e-10, '03': 1.759e-07, '02': 2.0295e-05}
        self.assertEqual(lines, exp_lines)
        self.assertEqual(significant, exp_significant)

    def test_get_data_enrichment_analysis_names(self):
        data = copy.deepcopy(ENRICH_DATA)
        data[LABEL] = LABEL_NAMES
        data, significant = get_data_enrichment_analysis(data, ENRICH_REF_AB, BINOMIAL_TEST)
        lines = data_to_lines(data)
        exp_lines = {(2, '02', 0, 'two', 25, 30, 0.5, 0.3, 0.3, 2.509702991379166),
                     (4, '04', 0, 'four', 1, 10, 0.02, 0.1, 0.1, -1.2341542222355069),
                     (9, '09', 2, 'nine', 1, 3, 0.02, 0.03, 0.03, 0.0),
                     (0, '00', '', 'r', 50, 100, 1, 1, 1, 0.0),
                     (6, '06', 1, 'six', nan, 5, nan, 0.05, 0.05, nan),
                     (7, '07', 1, 'seven', nan, 1, nan, 0.01, 0.01, nan),
                     (1, '01', 0, 'one', 5, 40, 0.1, 0.4, 0.4, -5.420266413988895),
                     (5, '05', 1, 'five', 5, 20, 0.1, 0.2, 0.2, -1.103304935668835),
                     (3, '03', 0, 'three', 20, 20, 0.4, 0.2, 0.2, 2.948803113091024),
                     (8, '08', 2, 'eight', 1, 1, 0.02, 0.01, 0.01, 0.4034095751193356)}
        exp_significant = {'01': 3.7996e-06, '03': 0.0011251149, '02': 0.0030924096}
        self.assertEqual(lines, exp_lines)
        self.assertEqual(significant, exp_significant)


# TOPOLOGY MANAGEMENT TESTS
# ==================================================================================================

MULTI_ROOT_DATA = {IDS: ['R', 'R-1', 'R-2', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'],
                   ONTO_ID: ['R', 'R-1', 'R-2', '00', '01', '02', '03', '04', '05', '06', '07',
                             '08', '09'],
                   PARENT: ['', 'R', 'R-1', 'R-2', '0', '0', '0', '0', '1', '1', '1', '2', '2'],
                   LABEL: ['R', 'R-1', 'R-2', '00', '01', '02', '03', '04', '05', '06', '07', '08',
                           '09'],
                   COUNT: [50, 50, 50, 50, 5, 25, 20, 1, 5, nan, nan, 1, 1],
                   REF_COUNT: [100, 100, 100, 100, 40, 30, 20, 10, 20, 5, 1, 1, 3],
                   PROP: [1, 1, 1, 1, 0.1, 0.5, 0.4, 0.02, 0.1, nan, nan, 0.02, 0.02],
                   REF_PROP: [1, 1, 1, 1, 0.4, 0.3, 0.2, 0.1, 0.2, 0.05, 0.01, 0.01, 0.03],
                   RELAT_PROP: [1000000, 1000000, 1000000, 1000000, 400000.0, 300000.0, 200000.0,
                                100000.0, 200000.0, 50000.0, 10000.0, 10000.0, 30000.0]}


class TestTopologyManagement(unittest.TestCase):

    @test_for(data_cut_root)
    def test_data_cut_root_uncut(self):
        data_b = copy.deepcopy(MULTI_ROOT_DATA)
        data = data_cut_root(data_b, ROOT_UNCUT)
        self.assertEqual(data, MULTI_ROOT_DATA)

    @test_for(data_cut_root)
    def test_data_cut_root_cut(self):
        data_b = copy.deepcopy(MULTI_ROOT_DATA)
        data = data_cut_root(data_b, ROOT_CUT)
        exp_data = {'ID': ['1', '2', '3', '4', '5', '6', '7', '8', '9'],
                    'Onto ID': ['01', '02', '03', '04', '05', '06', '07', '08', '09'],
                    'Parent': ['0', '0', '0', '0', '1', '1', '1', '2', '2'],
                    'Label': ['01', '02', '03', '04', '05', '06', '07', '08', '09'],
                    'Count': [5, 25, 20, 1, 5, nan, nan, 1, 1],
                    'Reference count': [40, 30, 20, 10, 20, 5, 1, 1, 3],
                    'Proportion': [0.1, 0.5, 0.4, 0.02, 0.1, nan, nan, 0.02, 0.02],
                    'Reference proportion': [0.4, 0.3, 0.2, 0.1, 0.2, 0.05, 0.01, 0.01, 0.03],
                    'Relative proportion': [400000.0, 300000.0, 200000.0, 100000.0, 200000.0,
                                            50000.0, 10000.0, 10000.0, 30000.0]}
        self.assertEqual(data, exp_data)

    @test_for(data_cut_root)
    def test_data_cut_root_total_cut(self):
        data_b = copy.deepcopy(MULTI_ROOT_DATA)
        data = data_cut_root(data_b, ROOT_TOTAL_CUT)
        exp_data = {'ID': ['1', '2', '3', '4', '5', '6', '7', '8', '9'],
                    'Onto ID': ['01', '02', '03', '04', '05', '06', '07', '08', '09'],
                    'Parent': ['', '', '', '', '1', '1', '1', '2', '2'],
                    'Label': ['01', '02', '03', '04', '05', '06', '07', '08', '09'],
                    'Count': [5, 25, 20, 1, 5, nan, nan, 1, 1],
                    'Reference count': [40, 30, 20, 10, 20, 5, 1, 1, 3],
                    'Proportion': [0.1, 0.5, 0.4, 0.02, 0.1, nan, nan, 0.02, 0.02],
                    'Reference proportion': [0.4, 0.3, 0.2, 0.1, 0.2, 0.05, 0.01, 0.01, 0.03],
                    'Relative proportion': [400000.0, 300000.0, 200000.0, 100000.0, 200000.0,
                                            50000.0, 10000.0, 10000.0, 30000.0]}
        self.assertEqual(data, exp_data)
