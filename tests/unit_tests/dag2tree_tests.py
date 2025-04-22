import unittest
import io

from functools import wraps
from ontosunburst.dag2tree import *

"""
Tests manually good file creation.
No automatic tests integrated.
"""

# ==================================================================================================
# GLOBAL
# ==================================================================================================

# --------------------------------------------------------------------------------------------------

ROOT = 'root'
CT_ONTO = {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf'], 'd': ['cde'], 'e': ['cde', 'eg'],
           'f': ['cf'], 'g': ['gh', 'eg'], 'h': ['gh'],
           'ab': [ROOT], 'cde': ['cdecf', 'cdeeg'], 'cf': ['cdecf'],
           'eg': [ROOT, 'cdeeg'], 'gh': [ROOT],
           'cdecf': [ROOT], 'cdeeg': ['cdeeg+'], 'cdeeg+': [ROOT]}

CT_AB = {ROOT: 6, 'cde': 3, 'cf': 3, 'cdecf': 3, 'cdeeg+': 3, 'cdeeg': 3, 'c': 3, 'ab': 3,
         'b': 2, 'a': 1}
CT_REF_AB = {ROOT: 36, 'cdeeg+': 19, 'cdeeg': 19, 'cdecf': 18, 'gh': 15, 'eg': 12, 'cde': 12,
             'cf': 9, 'h': 8, 'g': 7, 'f': 6, 'e': 5, 'd': 4, 'c': 3, 'ab': 3, 'b': 2, 'a': 1}

CT_LAB = {ROOT: 'Root', 'cdeeg+': 'CDEEG+', 'cdeeg': 'CDEEG', 'cdecf': 'CDECF', 'gh': 'GH',
          'eg': 'EG', 'cde': 'CDE', 'cf': 'CF', 'h': 'H', 'g': 'G', 'f': 'F', 'e': 'E', 'd': 'D',
          'c': 'C', 'ab': 'AB', 'b': 'B'}

W_PROP = {'1': 1.0, '2': 0.5, '3': 0.16666666666666666, '4': 0.3333333333333333, '5': nan,
          '6': nan, '7': nan, '8': nan, '9': nan, '10': nan, '11': 0.5, '12': 0.5, '13': 0.5,
          '14': nan, '15': nan, '16': 0.5, '17': 0.5, '18': nan, '19': 0.5, '20': 0.5, '21': 0.5,
          '22': 0.5, '23': nan, '24': nan, '25': nan, '26': nan, '27': nan}

W_REF_PROP = {'1': 1.0, '2': 0.08333333333333333, '3': 0.027777777777777776,
              '4': 0.05555555555555555, '5': 0.3333333333333333, '6': 0.1388888888888889,
              '7': 0.19444444444444445, '8': 0.4166666666666667, '9': 0.19444444444444445,
              '10': 0.2222222222222222, '11': 0.5, '12': 0.3333333333333333,
              '13': 0.08333333333333333, '14': 0.1111111111111111, '15': 0.1388888888888889,
              '16': 0.25, '17': 0.08333333333333333, '18': 0.16666666666666666,
              '19': 0.5277777777777778, '20': 0.5277777777777778, '21': 0.3333333333333333,
              '22': 0.08333333333333333, '23': 0.1111111111111111, '24': 0.1388888888888889,
              '25': 0.3333333333333333, '26': 0.1388888888888889, '27': 0.19444444444444445}

W_REL_PROP = {'1': 1000000, '2': 44776, '3': 14925, '4': 29850, '5': 179104, '6': 74626,
              '7': 104477, '8': 223880, '9': 104477, '10': 119402, '11': 268656, '12': 153517,
              '13': 38379, '14': 51172, '15': 63965, '16': 115138, '17': 38379, '18': 76758,
              '19': 283582, '20': 283582, '21': 141791, '22': 35447, '23': 47263, '24': 59079,
              '25': 141791, '26': 59079, '27': 82711}


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
        line = (dico[IDS][i], dico[ONTO_ID][i], dico[PARENT][i], dico[LABEL][i], dico[WEIGHT][i],
                dico[REF_WEIGHT][i])
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

    @test_for(get_set2_abundance)
    def test_get_sub_abundances_exists_diff(self):
        sub_abu = get_set2_abundance(CT_AB, 'cf')
        self.assertEqual(sub_abu, 3)

    @test_for(get_set2_abundance)
    def test_get_sub_abundances_exists_equ(self):
        sub_abu = get_set2_abundance(CT_AB, 'a')
        self.assertEqual(sub_abu, 1)

    @test_for(get_set2_abundance)
    def test_get_sub_abundances_not_exists(self):
        sub_abu = get_set2_abundance(CT_AB, 'eg')
        self.assertTrue(np.isnan(sub_abu))

    @test_for(TreeData.add_value)
    def test_add_value_data(self):
        data = TreeData()
        data.add_value(m_id='bjr', onto_id='Bjr_0', label='bonjour', count=2, ref_count=8,
                       parent='salutations')
        data.add_value(m_id='slt', onto_id='sl_1', label='salut', count=0.5, ref_count=2.3,
                       parent='salutations')
        wanted_data = {IDS: ['bjr', 'slt'],
                       ONTO_ID: ['Bjr_0', 'sl_1'],
                       PARENT: ['salutations', 'salutations'],
                       LABEL: ['bonjour', 'salut'],
                       WEIGHT: [2, 0.5],
                       REF_WEIGHT: [8, 2.3],
                       PROP: [nan, nan], REF_PROP: [nan, nan], RELAT_PROP: [nan, nan],
                       PVAL: [nan, nan]}
        self.assertEqual(data.get_data_dict(), wanted_data)

    @test_for(TreeData.dag_to_tree)
    def test_get_fig_parameters(self):
        data = TreeData()
        data.dag_to_tree(ref_abundance=CT_REF_AB, parent_dict=CT_ONTO,
                         root_item=ROOT, set_abundance=CT_AB, names=None)
        lines = set(data.get_col())
        w_lines = {('27', 'g', 'g', '25', nan, 7, nan, nan, nan, nan),
                   ('13', 'c', 'c', '12', 3, 3, nan, nan, nan, nan),
                   ('25', 'eg', 'eg', '20', nan, 12, nan, nan, nan, nan),
                   ('24', 'e', 'e', '21', nan, 5, nan, nan, nan, nan),
                   ('12', 'cde', 'cde', '11', 3, 12, nan, nan, nan, nan),
                   ('16', 'cf', 'cf', '11', 3, 9, nan, nan, nan, nan),
                   ('20', 'cdeeg', 'cdeeg', '19', 3, 19, nan, nan, nan, nan),
                   ('18', 'f', 'f', '16', nan, 6, nan, nan, nan, nan),
                   ('15', 'e', 'e', '12', nan, 5, nan, nan, nan, nan),
                   ('2', 'ab', 'ab', '1', 3, 3, nan, nan, nan, nan),
                   ('1', ROOT, ROOT, '', 6, 36, nan, nan, nan, nan),
                   ('8', 'gh', 'gh', '1', nan, 15, nan, nan, nan, nan),
                   ('26', 'e', 'e', '25', nan, 5, nan, nan, nan, nan),
                   ('19', 'cdeeg+', 'cdeeg+', '1', 3, 19, nan, nan, nan, nan),
                   ('17', 'c', 'c', '16', 3, 3, nan, nan, nan, nan),
                   ('23', 'd', 'd', '21', nan, 4, nan, nan, nan, nan),
                   ('11', 'cdecf', 'cdecf', '1', 3, 18, nan, nan, nan, nan),
                   ('7', 'g', 'g', '5', nan, 7, nan, nan, nan, nan),
                   ('9', 'g', 'g', '8', nan, 7, nan, nan, nan, nan),
                   ('4', 'b', 'b', '2', 2, 2, nan, nan, nan, nan),
                   ('5', 'eg', 'eg', '1', nan, 12, nan, nan, nan, nan),
                   ('3', 'a', 'a', '2', 1, 1, nan, nan, nan, nan),
                   ('22', 'c', 'c', '21', 3, 3, nan, nan, nan, nan),
                   ('21', 'cde', 'cde', '20', 3, 12, nan, nan, nan, nan),
                   ('14', 'd', 'd', '12', nan, 4, nan, nan, nan, nan),
                   ('10', 'h', 'h', '8', nan, 8, nan, nan, nan, nan),
                   ('6', 'e', 'e', '5', nan, 5, nan, nan, nan, nan)}
        self.assertEqual(lines, w_lines)

    @test_for(TreeData.dag_to_tree)
    def test_get_fig_parameters_names(self):
        data = TreeData()
        data.dag_to_tree(ref_abundance=CT_REF_AB, parent_dict=CT_ONTO,
                         root_item=ROOT, set_abundance=CT_AB, names=CT_LAB)
        lines = set(data.get_col())
        w_lines = {('9', 'g', 'G', '8', nan, 7, nan, nan, nan, nan),
                   ('19', 'cdeeg+', 'CDEEG+', '1', 3, 19, nan, nan, nan, nan),
                   ('14', 'd', 'D', '12', nan, 4, nan, nan, nan, nan),
                   ('23', 'd', 'D', '21', nan, 4, nan, nan, nan, nan),
                   ('6', 'e', 'E', '5', nan, 5, nan, nan, nan, nan),
                   ('25', 'eg', 'EG', '20', nan, 12, nan, nan, nan, nan),
                   ('16', 'cf', 'CF', '11', 3, 9, nan, nan, nan, nan),
                   ('10', 'h', 'H', '8', nan, 8, nan, nan, nan, nan),
                   ('8', 'gh', 'GH', '1', nan, 15, nan, nan, nan, nan),
                   ('1', ROOT, 'Root', '', 6, 36, nan, nan, nan, nan),
                   ('27', 'g', 'G', '25', nan, 7, nan, nan, nan, nan),
                   ('17', 'c', 'C', '16', 3, 3, nan, nan, nan, nan),
                   ('15', 'e', 'E', '12', nan, 5, nan, nan, nan, nan),
                   ('4', 'b', 'B', '2', 2, 2, nan, nan, nan, nan),
                   ('5', 'eg', 'EG', '1', nan, 12, nan, nan, nan, nan),
                   ('24', 'e', 'E', '21', nan, 5, nan, nan, nan, nan),
                   ('3', 'a', 'a', '2', 1, 1, nan, nan, nan, nan),
                   ('26', 'e', 'E', '25', nan, 5, nan, nan, nan, nan),
                   ('21', 'cde', 'CDE', '20', 3, 12, nan, nan, nan, nan),
                   ('20', 'cdeeg', 'CDEEG', '19', 3, 19, nan, nan, nan, nan),
                   ('13', 'c', 'C', '12', 3, 3, nan, nan, nan, nan),
                   ('12', 'cde', 'CDE', '11', 3, 12, nan, nan, nan, nan),
                   ('22', 'c', 'C', '21', 3, 3, nan, nan, nan, nan),
                   ('18', 'f', 'F', '16', nan, 6, nan, nan, nan, nan),
                   ('11', 'cdecf', 'CDECF', '1', 3, 18, nan, nan, nan, nan),
                   ('2', 'ab', 'AB', '1', 3, 3, nan, nan, nan, nan),
                   ('7', 'g', 'G', '5', nan, 7, nan, nan, nan, nan)}
        self.assertEqual(lines, w_lines)

    @test_for(TreeData.dag_to_tree)
    def test_get_fig_parameters_no_ref_base(self):
        data = TreeData()
        data.dag_to_tree(ref_abundance=CT_REF_AB, parent_dict=CT_ONTO,
                         root_item=ROOT, set_abundance=CT_AB, names=None,
                         ref_base=False)
        lines = set(data.get_col())
        w_lines = {('4', 'b', 'b', '2', 2, 2, nan, nan, nan, nan),
                   ('7', 'c', 'c', '6', 3, 3, nan, nan, nan, nan),
                   ('13', 'c', 'c', '12', 3, 3, nan, nan, nan, nan),
                   ('8', 'cf', 'cf', '5', 3, 9, nan, nan, nan, nan),
                   ('3', 'a', 'a', '2', 1, 1, nan, nan, nan, nan),
                   ('2', 'ab', 'ab', '1', 3, 3, nan, nan, nan, nan),
                   ('9', 'c', 'c', '8', 3, 3, nan, nan, nan, nan),
                   ('10', 'cdeeg+', 'cdeeg+', '1', 3, 19, nan, nan, nan, nan),
                   ('12', 'cde', 'cde', '11', 3, 12, nan, nan, nan, nan),
                   ('5', 'cdecf', 'cdecf', '1', 3, 18, nan, nan, nan, nan),
                   ('6', 'cde', 'cde', '5', 3, 12, nan, nan, nan, nan),
                   ('1', ROOT, ROOT, '', 6, 36, nan, nan, nan, nan),
                   ('11', 'cdeeg', 'cdeeg', '10', 3, 19, nan, nan, nan, nan)}
        self.assertEqual(lines, w_lines)

    @test_for(TreeData.dag_to_tree)
    def test_get_fig_parameters_no_ref(self):
        data = TreeData()
        data.dag_to_tree(ref_abundance=CT_AB, parent_dict=CT_ONTO,
                         root_item=ROOT, set_abundance=CT_AB, names=None,
                         ref_base=False)
        lines = set(data.get_col())
        w_lines = {('1', ROOT, ROOT, '', 6, 6, nan, nan, nan, nan),
                   ('8', 'cf', 'cf', '5', 3, 3, nan, nan, nan, nan),
                   ('13', 'c', 'c', '12', 3, 3, nan, nan, nan, nan),
                   ('9', 'c', 'c', '8', 3, 3, nan, nan, nan, nan),
                   ('11', 'cdeeg', 'cdeeg', '10', 3, 3, nan, nan, nan, nan),
                   ('7', 'c', 'c', '6', 3, 3, nan, nan, nan, nan),
                   ('3', 'a', 'a', '2', 1, 1, nan, nan, nan, nan),
                   ('12', 'cde', 'cde', '11', 3, 3, nan, nan, nan, nan),
                   ('2', 'ab', 'ab', '1', 3, 3, nan, nan, nan, nan),
                   ('10', 'cdeeg+', 'cdeeg+', '1', 3, 3, nan, nan, nan, nan),
                   ('6', 'cde', 'cde', '5', 3, 3, nan, nan, nan, nan),
                   ('5', 'cdecf', 'cdecf', '1', 3, 3, nan, nan, nan, nan),
                   ('4', 'b', 'b', '2', 2, 2, nan, nan, nan, nan)}
        self.assertEqual(lines, w_lines)


class TestAddProportionDataTable(unittest.TestCase):

    @test_for(TreeData.calculate_proportions)
    def test_get_data_proportion_no_relative(self):
        data = TreeData()
        data.dag_to_tree(ref_abundance=CT_REF_AB, parent_dict=CT_ONTO,
                         root_item=ROOT, set_abundance=CT_AB, names=CT_LAB)
        data.calculate_proportions(True)
        for i in range(data.len):
            if np.isnan(data.prop[i]):
                self.assertTrue(np.isnan(W_PROP[data.ids[i]]))
            else:
                self.assertEqual(data.prop[i], W_PROP[data.ids[i]])

    @test_for(TreeData.calculate_proportions)
    def test_get_data_proportion_no_relative_ref(self):
        data = TreeData()
        data.dag_to_tree(ref_abundance=CT_REF_AB, parent_dict=CT_ONTO,
                         root_item=ROOT, set_abundance=CT_AB, names=CT_LAB)
        data.calculate_proportions(True)
        for i in range(data.len):
            self.assertEqual(data.ref_prop[i], W_REF_PROP[data.ids[i]])

    @test_for(TreeData.calculate_proportions)
    def test_get_data_proportion_relative(self):
        data = TreeData()
        data.dag_to_tree(ref_abundance=CT_REF_AB, parent_dict=CT_ONTO,
                         root_item=ROOT, set_abundance=CT_AB, names=CT_LAB)
        data.calculate_proportions(True)
        for i in range(data.len):
            self.assertEqual(data.relative_prop[i], W_REL_PROP[data.ids[i]])


# ENRICHMENT TESTS
# ==================================================================================================

ENRICH_AB = {'00': 50, '01': 5, '02': 25, '03': 20, '04': 1, '05': 5, '06': nan, '07': nan,
             '08': 1, '09': 1}
ENRICH_REF_AB = {'00': 100, '01': 40, '02': 30, '03': 20, '04': 10, '05': 20, '06': 5, '07': 1,
                 '08': 1, '09': 3}
E_LABElS = {'00': '0', '01': '1', '02': '2', '03': '3', '04': '4',
            '05': '5', '06': '6', '07': '7', '08': '8', '09': '9'}
E_ONTO = {'01': ['00'], '02': ['00'], '03': ['00'], '04': ['00'], '05': ['01'],
          '06': ['01'], '07': ['01'], '08': ['02'], '09': ['02']}


# Expected :
# Over : 2, 3 | Under : 1, 4, 5 | No diff : 0, 8, 9 | Nan : 6, 7


class TestEnrichmentAnalysis(unittest.TestCase):

    @test_for(TreeData.make_enrichment_analysis)
    def test_get_data_enrichment_analysis_single_value(self):
        data = TreeData()
        data.dag_to_tree(ENRICH_AB, ENRICH_REF_AB, E_ONTO, '00', E_LABElS)
        data.calculate_proportions(True)
        data.make_enrichment_analysis(BINOMIAL_TEST)
        p_value_1 = [data.p_val[i] for i in range(data.len) if data.onto_ids[i] == '01'][0]
        M = 100
        N = 50
        m = 40
        n = 5
        exp_p_value_1 = stats.binomtest(n, N, m / M, alternative='two-sided').pvalue
        exp_p_value_1 = np.log10(exp_p_value_1)
        self.assertEqual(p_value_1, exp_p_value_1)

    @test_for(TreeData.make_enrichment_analysis)
    def test_get_data_enrichment_analysis_binomial(self):
        data = TreeData()
        data.dag_to_tree(ENRICH_AB, ENRICH_REF_AB, E_ONTO, '00', E_LABElS)
        data.calculate_proportions(True)
        significant = data.make_enrichment_analysis(BINOMIAL_TEST)
        lines = set(data.get_col())
        print(lines)
        exp_significant = {'01': 3.799562441228011e-06, '03': 0.001125114927936431,
                           '02': 0.003092409570144631}
        exp_lines = {('8', '09', '9', '6', 1, 3, 0.02, 0.03, 30000, 0.0),
                     ('9', '03', '3', '1', 20, 20, 0.4, 0.2, 200000, 2.948803113091024),
                     ('5', '07', '7', '2', nan, 1, nan, 0.01, 10000, nan),
                     ('2', '01', '1', '1', 5, 40, 0.1, 0.4, 400000, -5.420266413988895),
                     ('6', '02', '2', '1', 25, 30, 0.5, 0.3, 300000, 2.509702991379166),
                     ('4', '06', '6', '2', nan, 5, nan, 0.05, 50000, nan),
                     ('7', '08', '8', '6', 1, 1, 0.02, 0.01, 10000, 0.4034095751193356),
                     ('10', '04', '4', '1', 1, 10, 0.02, 0.1, 100000, -1.2341542222355069),
                     ('1', '00', '0', '', 50, 100, 1.0, 1.0, 1000000, 0.0),
                     ('3', '05', '5', '2', 5, 20, 0.1, 0.2, 200000, -1.103304935668835)}

        for line in lines:
            line = tuple([nan if type(x) != str and np.isnan(x) else x for x in line])
            self.assertIn(line, exp_lines)
        self.assertEqual(len(lines), len(exp_lines))
        self.assertEqual(significant, exp_significant)

    @test_for(TreeData.make_enrichment_analysis)
    def test_get_data_enrichment_analysis_hypergeometric(self):
        data = TreeData()
        data.dag_to_tree(ENRICH_AB, ENRICH_REF_AB, E_ONTO, '00', E_LABElS)
        data.calculate_proportions(True)
        significant = data.make_enrichment_analysis(HYPERGEO_TEST)
        lines = set(data.get_col())
        exp_lines = {('6', '02', '2', '1', 25, 30, 0.5, 0.3, 300000, 4.692610428021241),
                     ('2', '01', '1', '1', 5, 40, 0.1, 0.4, 400000, -9.138873998573988),
                     ('1', '00', '0', '', 50, 100, 1.0, 1.0, 1000000, 0.3010299956639812),
                     ('10', '04', '4', '1', 1, 10, 0.02, 0.1, 100000, -1.8051946563380086),
                     ('3', '05', '5', '2', 5, 20, 0.1, 0.2, 200000, -1.6413993451973743),
                     ('9', '03', '3', '1', 20, 20, 0.4, 0.2, 200000, 6.754831139005899),
                     ('4', '06', '6', '2', nan, 5, nan, 0.05, 50000, nan),
                     ('8', '09', '9', '6', 1, 3, 0.02, 0.03, 30000, -1.4464911998299308e-16),
                     ('5', '07', '7', '2', nan, 1, nan, 0.01, 10000, nan),
                     ('7', '08', '8', '6', 1, 1, 0.02, 0.01, 10000, -0.0)}
        exp_significant = {'01': 7.263166523971598e-10, '03': 1.7586072571039978e-07,
                           '02': 2.0295024128400847e-05}
        for line in lines:
            line = tuple([nan if type(x) != str and np.isnan(x) else x for x in line])
            self.assertIn(line, exp_lines)
        self.assertEqual(len(lines), len(exp_lines))
        self.assertEqual(significant, exp_significant)

    @test_for(TreeData.make_enrichment_analysis)
    def test_get_data_enrichment_analysis_scores(self):
        scores = {'00': 0.05, '01': 0.2, '02': 0.0004, '03': 0.5, '04': 0.000008, '05': 0.9,
                  '06': 0.01, '07': nan, '08': nan, '09': 0.000023}
        data = TreeData()
        data.dag_to_tree(ENRICH_AB, ENRICH_REF_AB, E_ONTO, '00', E_LABElS)
        data.calculate_proportions(True)
        significant = data.make_enrichment_analysis(HYPERGEO_TEST, scores)
        lines = set(data.get_col())
        exp_lines = {('9', '03', '3', '1', 20, 20, 0.4, 0.2, 200000, 0.3010299956639812),
                     ('8', '09', '9', '6', 1, 3, 0.02, 0.03, 30000, 4.638272163982407),
                     ('2', '01', '1', '1', 5, 40, 0.1, 0.4, 400000, 0.6989700043360187),
                     ('3', '05', '5', '2', 5, 20, 0.1, 0.2, 200000, 0.045757490560675115),
                     ('4', '06', '6', '2', nan, 5, nan, 0.05, 50000, 2.0),
                     ('6', '02', '2', '1', 25, 30, 0.5, 0.3, 300000, 3.3979400086720375),
                     ('5', '07', '7', '2', nan, 1, nan, 0.01, 10000, nan),
                     ('7', '08', '8', '6', 1, 1, 0.02, 0.01, 10000, nan),
                     ('1', '00', '0', '', 50, 100, 1.0, 1.0, 1000000, 1.3010299956639813),
                     ('10', '04', '4', '1', 1, 10, 0.02, 0.1, 100000, 5.096910013008056)}
        exp_significant = {'04': 8e-06, '09': 2.3e-05, '02': 0.0004}
        for line in lines:
            line = tuple([nan if type(x) != str and np.isnan(x) else x for x in line])
            self.assertIn(line, exp_lines)
        self.assertEqual(len(lines), len(exp_lines))
        self.assertEqual(significant, exp_significant)


# TOPOLOGY MANAGEMENT TESTS
# ==================================================================================================

ROOT_AB = {'R': 50, 'R-1': 50, 'R-2': 50, '00': 50, '01': 5, '02': 25, '03': 20, '04': 1, '05': 5,
           '06': nan, '07': nan, '08': 1, '09': 1}
ROOT_REF_AB = {'R': 100, 'R-1': 100, 'R-2': 100, '00': 100, '01': 40, '02': 30, '03': 20, '04': 10,
               '05': 20, '06': 5, '07': 1, '08': 1, '09': 3}
ROOT_LABElS = {'R': 'r', '00': '0', '01': '1', '02': '2', '03': '3',
               '04': '4', '05': '5', '06': '6', '07': '7', '08': '8', '09': '9'}
ROOT_LABElS_C = {'R': 'r', '00': '6', '01': '1', '02': '2', '03': '3',
                 '04': '4', '05': '5', '06': '6', '07': '7', '08': '8', '09': '9'}
ROOT_ONTO = {'01': ['00'], '02': ['00'], '03': ['00'], '04': ['00'], '05': ['01'], '00': ['R-2'],
             '06': ['01'], '07': ['01'], '08': ['02'], '09': ['02'], 'R-1': ['R'], 'R-2': ['R-1']}

PATH_ONTO = {'a': ['ab', 'cdeeg++++'], 'b': ['ab'], 'c': ['cde', 'cf'], 'd': ['cde'],
             'e': ['cde', 'eg'], 'f': ['cf'], 'g': ['gh', 'eg'], 'h': ['gh'],
             'ab': [ROOT], 'cde': ['cde+'], 'cde+': ['cde++'], 'cde++': ['cde+++'],
             'cde+++': ['cdecf', 'cdeeg'], 'cf': ['cdecf'], 'eg': ['cdeeg', ROOT],
             'gh': [ROOT], 'cdecf': [ROOT], 'cdeeg': ['cdeeg+'],
             'cdeeg+': ['cdeeg++'], 'cdeeg++': ['cdeeg+++'], 'cdeeg+++': ['cdeeg++++'],
             'cdeeg++++': [ROOT]}
PATH_AB = {ROOT: 6, 'cde': 3, 'cde+': 3, 'cde++': 3, 'cde+++': 3, 'cf': 3, 'cdecf': 3,
           'cdeeg++++': 3, 'cdeeg+++': 3, 'cdeeg++': 3, 'cdeeg+': 3, 'cdeeg': 3, 'c': 3, 'ab': 3,
           'b': 2, 'a': 1}
PATH_REF_AB = {ROOT: 36, 'cdeeg++++': 19, 'cdeeg+++': 19, 'cdeeg++': 19, 'cdeeg+': 19,
               'cdeeg': 19, 'cdecf': 18, 'gh': 15, 'eg': 12, 'cde': 12, 'cde+': 12, 'cde++': 12,
               'cde+++': 12, 'cf': 9, 'h': 8, 'g': 7, 'f': 6, 'e': 5, 'd': 4, 'c': 3, 'ab': 3,
               'b': 2, 'a': 1}
PATH_LAB = {ROOT: 'Root', 'cdeeg+': 'CDEEG+', 'cdeeg': 'CDEEG', 'cdecf': 'CDECF', 'gh': 'GH',
            'eg': 'EG', 'cde': 'CDE', 'cf': 'CF', 'h': 'H', 'g': 'G', 'f': 'F', 'e': 'E', 'd': 'D',
            'c': 'C', 'ab': 'AB', 'b': 'B'}


class TestTopologyManagement(unittest.TestCase):

    @test_for(TreeData.cut_root)
    def test_data_cut_root_uncut(self):
        data = TreeData()
        data.dag_to_tree(ROOT_AB, ROOT_REF_AB, ROOT_ONTO, 'R', ROOT_LABElS)
        data.calculate_proportions(True)
        exp_d = data.get_data_dict()
        data.cut_root(ROOT_UNCUT)
        self.assertEqual(data.get_data_dict(), exp_d)

    @test_for(TreeData.cut_root)
    def test_data_cut_root_cut(self):
        data = TreeData()
        data.dag_to_tree(ROOT_AB, ROOT_REF_AB, ROOT_ONTO, 'R', ROOT_LABElS)
        data.calculate_proportions(True)
        data.cut_root(ROOT_CUT)
        lines = set(data.get_col())
        exp_lines = {('12', '03', '3', '0', 20, 20, 0.4, 0.2, 200000, nan),
                     ('13', '04', '4', '0', 1, 10, 0.02, 0.1, 100000, nan),
                     ('8', '07', '7', '5', nan, 1, nan, 0.01, 10000, nan),
                     ('11', '09', '9', '9', 1, 3, 0.02, 0.03, 30000, nan),
                     ('9', '02', '2', '0', 25, 30, 0.5, 0.3, 300000, nan),
                     ('5', '01', '1', '0', 5, 40, 0.1, 0.4, 400000, nan),
                     ('6', '05', '5', '5', 5, 20, 0.1, 0.2, 200000, nan),
                     ('7', '06', '6', '5', nan, 5, nan, 0.05, 50000, nan),
                     ('10', '08', '8', '9', 1, 1, 0.02, 0.01, 10000, nan)}
        self.assertEqual(data.len, 9)
        self.assertEqual(len(lines), len(exp_lines))
        for line in lines:
            line = tuple([nan if type(x) != str and np.isnan(x) else x for x in line])
            self.assertIn(line, exp_lines)

    @test_for(TreeData.cut_root)
    def test_data_cut_root_cut_id_conflict(self):
        data = TreeData()
        data.dag_to_tree(ROOT_AB, ROOT_REF_AB, ROOT_ONTO, 'R', ROOT_LABElS_C)
        data.calculate_proportions(True)
        data.cut_root(ROOT_CUT)
        lines = set(data.get_col())
        exp_lines = {('9', '02', '2', '6_', 25, 30, 0.5, 0.3, 300000, nan),
                     ('12', '03', '3', '6_', 20, 20, 0.4, 0.2, 200000, nan),
                     ('5', '01', '1', '6_', 5, 40, 0.1, 0.4, 400000, nan),
                     ('11', '09', '9', '9', 1, 3, 0.02, 0.03, 30000, nan),
                     ('7', '06', '6', '5', nan, 5, nan, 0.05, 50000, nan),
                     ('13', '04', '4', '6_', 1, 10, 0.02, 0.1, 100000, nan),
                     ('6', '05', '5', '5', 5, 20, 0.1, 0.2, 200000, nan),
                     ('10', '08', '8', '9', 1, 1, 0.02, 0.01, 10000, nan),
                     ('8', '07', '7', '5', nan, 1, nan, 0.01, 10000, nan)}

        self.assertEqual(data.len, 9)
        self.assertEqual(len(lines), len(exp_lines))
        for line in lines:
            line = tuple([nan if type(x) != str and np.isnan(x) else x for x in line])
            self.assertIn(line, exp_lines)

    @test_for(TreeData.cut_root)
    def test_data_cut_root_total_cut(self):
        data = TreeData()
        data.dag_to_tree(ROOT_AB, ROOT_REF_AB, ROOT_ONTO, 'R', ROOT_LABElS)
        data.calculate_proportions(True)
        data.cut_root(ROOT_TOTAL_CUT)
        lines = set(data.get_col())
        exp_lines = {('9', '02', '2', '', 25, 30, 0.5, 0.3, 300000, nan),
                     ('6', '05', '5', '5', 5, 20, 0.1, 0.2, 200000, nan),
                     ('10', '08', '8', '9', 1, 1, 0.02, 0.01, 10000, nan),
                     ('13', '04', '4', '', 1, 10, 0.02, 0.1, 100000, nan),
                     ('5', '01', '1', '', 5, 40, 0.1, 0.4, 400000, nan),
                     ('11', '09', '9', '9', 1, 3, 0.02, 0.03, 30000, nan),
                     ('7', '06', '6', '5', nan, 5, nan, 0.05, 50000, nan),
                     ('12', '03', '3', '', 20, 20, 0.4, 0.2, 200000, nan),
                     ('8', '07', '7', '5', nan, 1, nan, 0.01, 10000, nan)}
        self.assertEqual(data.len, 9)
        self.assertEqual(len(lines), len(exp_lines))
        for line in lines:
            line = tuple([nan if type(x) != str and np.isnan(x) else x for x in line])
            self.assertIn(line, exp_lines)

    @test_for(TreeData.cut_nested_path)
    def test_cut_path_uncut(self):
        data = TreeData()
        data.dag_to_tree(PATH_AB, PATH_REF_AB, PATH_ONTO, ROOT, PATH_LAB)
        data.calculate_proportions(True)
        data.cut_nested_path(PATH_UNCUT, False)
        exp_l = {('27', 'cdeeg', 'CDEEG', '26', 3, 19, 0.5, 0.5277777777777778, 269402, nan),
                 ('21', 'f', 'F', '19', nan, 6, nan, 0.16666666666666666, 76758, nan),
                 ('4', 'b', 'B', '2', 2, 2, 0.3333333333333333, 0.05555555555555555, 29850, nan),
                 ('24', 'cdeeg+++', 'cdeeg+++', '22', 3, 19, 0.5, 0.5277777777777778, 269402, nan),
                 ('6', 'e', 'E', '5', nan, 5, nan, 0.1388888888888889, 74626, nan),
                 ('9', 'g', 'G', '8', nan, 7, nan, 0.19444444444444445, 104477, nan),
                 ('22', 'cdeeg++++', 'cdeeg++++', '1', 3, 19, 0.5, 0.5277777777777778, 283582, nan),
                 ('1', ROOT, 'Root', '', 6, 36, 1.0, 1.0, 1000000, nan),
                 ('20', 'c', 'C', '19', 3, 3, 0.5, 0.08333333333333333, 38379, nan),
                 ('8', 'gh', 'GH', '1', nan, 15, nan, 0.4166666666666667, 223880, nan),
                 (
                 '23', 'a', 'a', '22', 1, 1, 0.16666666666666666, 0.027777777777777776, 14179, nan),
                 ('16', 'c', 'C', '15', 3, 3, 0.5, 0.08333333333333333, 38379, nan),
                 ('36', 'e', 'E', '35', nan, 5, nan, 0.1388888888888889, 56125, nan),
                 ('14', 'cde+', 'cde+', '13', 3, 12, 0.5, 0.3333333333333333, 153517, nan),
                 ('15', 'cde', 'CDE', '14', 3, 12, 0.5, 0.3333333333333333, 153517, nan),
                 ('32', 'c', 'C', '31', 3, 3, 0.5, 0.08333333333333333, 33675, nan),
                 ('34', 'e', 'E', '31', nan, 5, nan, 0.1388888888888889, 56125, nan),
                 ('2', 'ab', 'AB', '1', 3, 3, 0.5, 0.08333333333333333, 44776, nan),
                 ('7', 'g', 'G', '5', nan, 7, nan, 0.19444444444444445, 104477, nan),
                 ('30', 'cde+', 'cde+', '29', 3, 12, 0.5, 0.3333333333333333, 134701, nan),
                 ('17', 'd', 'D', '15', nan, 4, nan, 0.1111111111111111, 51172, nan),
                 ('18', 'e', 'E', '15', nan, 5, nan, 0.1388888888888889, 63965, nan),
                 ('19', 'cf', 'CF', '11', 3, 9, 0.5, 0.25, 115138, nan),
                 ('33', 'd', 'D', '31', nan, 4, nan, 0.1111111111111111, 44900, nan),
                 ('37', 'g', 'G', '35', nan, 7, nan, 0.19444444444444445, 78575, nan),
                 ('13', 'cde++', 'cde++', '12', 3, 12, 0.5, 0.3333333333333333, 153517, nan),
                 ('25', 'cdeeg++', 'cdeeg++', '24', 3, 19, 0.5, 0.5277777777777778, 269402, nan),
                 ('3', 'a', 'a', '2', 1, 1, 0.16666666666666666, 0.027777777777777776, 14925, nan),
                 ('11', 'cdecf', 'CDECF', '1', 3, 18, 0.5, 0.5, 268656, nan),
                 ('10', 'h', 'H', '8', nan, 8, nan, 0.2222222222222222, 119402, nan),
                 ('28', 'cde+++', 'cde+++', '27', 3, 12, 0.5, 0.3333333333333333, 134701, nan),
                 ('31', 'cde', 'CDE', '30', 3, 12, 0.5, 0.3333333333333333, 134701, nan),
                 ('29', 'cde++', 'cde++', '28', 3, 12, 0.5, 0.3333333333333333, 134701, nan),
                 ('35', 'eg', 'EG', '27', nan, 12, nan, 0.3333333333333333, 134701, nan),
                 ('26', 'cdeeg+', 'CDEEG+', '25', 3, 19, 0.5, 0.5277777777777778, 269402, nan),
                 ('5', 'eg', 'EG', '1', nan, 12, nan, 0.3333333333333333, 179104, nan),
                 ('12', 'cde+++', 'cde+++', '11', 3, 12, 0.5, 0.3333333333333333, 153517, nan)}
        lines = set(data.get_col())
        self.assertEqual(len(lines), len(exp_l))
        for line in lines:
            line = tuple([nan if type(x) != str and np.isnan(x) else x for x in line])
            self.assertIn(line, exp_l)
        # from ontosunburst.sunburst_fig import generate_sunburst_fig
        # generate_sunburst_fig(data, 'test', bg_color='black', font_color='white')

    @test_for(TreeData.cut_nested_path)
    def test_cut_path_cut_deeper(self):
        data = TreeData()
        data.dag_to_tree(PATH_AB, PATH_REF_AB, PATH_ONTO, ROOT, PATH_LAB)
        data.calculate_proportions(True)
        data.cut_nested_path(PATH_DEEPER, False)
        # from ontosunburst.sunburst_fig import generate_sunburst_fig
        # generate_sunburst_fig(data, 'test')
        exp_l = {('34', 'e', 'E', '31', nan, 5, nan, 0.1388888888888889, 56125, nan),
                 ('1', ROOT, 'Root', '', 6, 36, 1.0, 1.0, 1000000, nan),
                 ('8', 'gh', 'GH', '1', nan, 15, nan, 0.4166666666666667, 223880, nan),
                 ('19', 'cf', 'CF', '11', 3, 9, 0.5, 0.25, 115138, nan),
                 ('27', 'cdeeg', '... CDEEG', '22', 3, 19, 0.5, 0.5277777777777778, 269402, nan),
                 ('2', 'ab', 'AB', '1', 3, 3, 0.5, 0.08333333333333333, 44776, nan),
                 ('32', 'c', 'C', '31', 3, 3, 0.5, 0.08333333333333333, 33675, nan),
                 ('11', 'cdecf', 'CDECF', '1', 3, 18, 0.5, 0.5, 268656, nan),
                 ('22', 'cdeeg++++', 'cdeeg++++', '1', 3, 19, 0.5, 0.5277777777777778, 283582, nan),
                 ('15', 'cde', '... CDE', '11', 3, 12, 0.5, 0.3333333333333333, 153517, nan),
                 ('5', 'eg', 'EG', '1', nan, 12, nan, 0.3333333333333333, 179104, nan),
                 ('37', 'g', 'G', '35', nan, 7, nan, 0.19444444444444445, 78575, nan),
                 ('16', 'c', 'C', '15', 3, 3, 0.5, 0.08333333333333333, 38379, nan),
                 ('33', 'd', 'D', '31', nan, 4, nan, 0.1111111111111111, 44900, nan),
                 ('36', 'e', 'E', '35', nan, 5, nan, 0.1388888888888889, 56125, nan),
                 ('35', 'eg', 'EG', '27', nan, 12, nan, 0.3333333333333333, 134701, nan),
                 ('7', 'g', 'G', '5', nan, 7, nan, 0.19444444444444445, 104477, nan),
                 ('9', 'g', 'G', '8', nan, 7, nan, 0.19444444444444445, 104477, nan),
                 ('4', 'b', 'B', '2', 2, 2, 0.3333333333333333, 0.05555555555555555, 29850, nan),
                 ('3', 'a', 'a', '2', 1, 1, 0.16666666666666666, 0.027777777777777776, 14925, nan),
                 (
                 '23', 'a', 'a', '22', 1, 1, 0.16666666666666666, 0.027777777777777776, 14179, nan),
                 ('6', 'e', 'E', '5', nan, 5, nan, 0.1388888888888889, 74626, nan),
                 ('20', 'c', 'C', '19', 3, 3, 0.5, 0.08333333333333333, 38379, nan),
                 ('10', 'h', 'H', '8', nan, 8, nan, 0.2222222222222222, 119402, nan),
                 ('31', 'cde', '... CDE', '27', 3, 12, 0.5, 0.3333333333333333, 134701, nan),
                 ('18', 'e', 'E', '15', nan, 5, nan, 0.1388888888888889, 63965, nan),
                 ('17', 'd', 'D', '15', nan, 4, nan, 0.1111111111111111, 51172, nan),
                 ('21', 'f', 'F', '19', nan, 6, nan, 0.16666666666666666, 76758, nan)}
        lines = set(data.get_col())
        self.assertEqual(len(lines), len(exp_l))
        for line in lines:
            line = tuple([nan if type(x) != str and np.isnan(x) else x for x in line])
            self.assertIn(line, exp_l)

    @test_for(TreeData.cut_nested_path)
    def test_cut_path_cut_higher(self):
        data = TreeData()
        data.dag_to_tree(PATH_AB, PATH_REF_AB, PATH_ONTO, ROOT, PATH_LAB)
        data.calculate_proportions(True)
        data.cut_root(ROOT_CUT)
        data.cut_nested_path(PATH_HIGHER, False)
        exp_l = {('10', 'h', 'H', '8', nan, 8, nan, 0.2222222222222222, 119402, nan),
                 ('11', 'cdecf', 'CDECF', 'Root', 3, 18, 0.5, 0.5, 268656, nan),
                 ('17', 'd', 'D', '12', nan, 4, nan, 0.1111111111111111, 51172, nan),
                 ('12', 'cde+++', 'cde+++ ...', '11', 3, 12, 0.5, 0.3333333333333333, 153517, nan),
                 ('28', 'cde+++', 'cde+++ ...', '24', 3, 12, 0.5, 0.3333333333333333, 134701, nan),
                 ('21', 'f', 'F', '19', nan, 6, nan, 0.16666666666666666, 76758, nan),
                 ('6', 'e', 'E', '5', nan, 5, nan, 0.1388888888888889, 74626, nan),
                 ('23', 'a', 'a', '22', 1, 1, 0.16666666666666666, 0.027777777777777776, 14179, nan),
                 ('20', 'c', 'C', '19', 3, 3, 0.5, 0.08333333333333333, 38379, nan),
                 ('32', 'c', 'C', '28', 3, 3, 0.5, 0.08333333333333333, 33675, nan),
                 ('4', 'b', 'B', '2', 2, 2, 0.3333333333333333, 0.05555555555555555, 29850, nan),
                 ('8', 'gh', 'GH', 'Root', nan, 15, nan, 0.4166666666666667, 223880, nan),
                 ('2', 'ab', 'AB', 'Root', 3, 3, 0.5, 0.08333333333333333, 44776, nan),
                 ('33', 'd', 'D', '28', nan, 4, nan, 0.1111111111111111, 44900, nan),
                 ('7', 'g', 'G', '5', nan, 7, nan, 0.19444444444444445, 104477, nan),
                 ('19', 'cf', 'CF', '11', 3, 9, 0.5, 0.25, 115138, nan),
                 ('22', 'cdeeg++++', 'cdeeg++++', 'Root', 3, 19, 0.5, 0.5277777777777778, 283582, nan),
                 ('16', 'c', 'C', '12', 3, 3, 0.5, 0.08333333333333333, 38379, nan),
                 ('34', 'e', 'E', '28', nan, 5, nan, 0.1388888888888889, 56125, nan),
                 ('3', 'a', 'a', '2', 1, 1, 0.16666666666666666, 0.027777777777777776, 14925, nan),
                 ('18', 'e', 'E', '12', nan, 5, nan, 0.1388888888888889, 63965, nan),
                 ('24', 'cdeeg+++', 'cdeeg+++ ...', '22', 3, 19, 0.5, 0.5277777777777778, 269402, nan),
                 ('35', 'eg', 'EG', '24', nan, 12, nan, 0.3333333333333333, 134701, nan),
                 ('36', 'e', 'E', '35', nan, 5, nan, 0.1388888888888889, 56125, nan),
                 ('37', 'g', 'G', '35', nan, 7, nan, 0.19444444444444445, 78575, nan),
                 ('5', 'eg', 'EG', 'Root', nan, 12, nan, 0.3333333333333333, 179104, nan),
                 ('9', 'g', 'G', '8', nan, 7, nan, 0.19444444444444445, 104477, nan)}
        lines = set(data.get_col())
        # from ontosunburst.sunburst_fig import generate_sunburst_fig
        # generate_sunburst_fig(data, 'test')
        self.assertEqual(len(lines), len(exp_l))
        for line in lines:
            line = tuple([nan if type(x) != str and np.isnan(x) else x for x in line])
            self.assertIn(line, exp_l)

    @test_for(TreeData.cut_nested_path)
    def test_cut_path_cut_bound(self):
        data = TreeData()
        data.dag_to_tree(PATH_AB, PATH_REF_AB, PATH_ONTO, ROOT, PATH_LAB)
        data.calculate_proportions(True)
        data.cut_nested_path(PATH_BOUND, False)
        exp_l = {('16', 'c', 'C', '15', 3, 3, 0.5, 0.08333333333333333, 38379, nan),
                 ('20', 'c', 'C', '19', 3, 3, 0.5, 0.08333333333333333, 38379, nan),
                 ('8', 'gh', 'GH', '1', nan, 15, nan, 0.4166666666666667, 223880, nan),
                 (
                 '23', 'a', 'a', '22', 1, 1, 0.16666666666666666, 0.027777777777777776, 14179, nan),
                 ('18', 'e', 'E', '15', nan, 5, nan, 0.1388888888888889, 63965, nan),
                 ('27', 'cdeeg', '... CDEEG', '24', 3, 19, 0.5, 0.5277777777777778, 269402, nan),
                 ('34', 'e', 'E', '31', nan, 5, nan, 0.1388888888888889, 56125, nan),
                 ('15', 'cde', '... CDE', '12', 3, 12, 0.5, 0.3333333333333333, 153517, nan),
                 ('6', 'e', 'E', '5', nan, 5, nan, 0.1388888888888889, 74626, nan),
                 ('32', 'c', 'C', '31', 3, 3, 0.5, 0.08333333333333333, 33675, nan),
                 ('21', 'f', 'F', '19', nan, 6, nan, 0.16666666666666666, 76758, nan),
                 ('7', 'g', 'G', '5', nan, 7, nan, 0.19444444444444445, 104477, nan),
                 ('36', 'e', 'E', '35', nan, 5, nan, 0.1388888888888889, 56125, nan),
                 ('11', 'cdecf', 'CDECF', '1', 3, 18, 0.5, 0.5, 268656, nan),
                 ('22', 'cdeeg++++', 'cdeeg++++', '1', 3, 19, 0.5, 0.5277777777777778, 283582, nan),
                 ('10', 'h', 'H', '8', nan, 8, nan, 0.2222222222222222, 119402, nan),
                 ('9', 'g', 'G', '8', nan, 7, nan, 0.19444444444444445, 104477, nan),
                 ('3', 'a', 'a', '2', 1, 1, 0.16666666666666666, 0.027777777777777776, 14925, nan),
                 ('12', 'cde+++', 'cde+++ ...', '11', 3, 12, 0.5, 0.3333333333333333, 153517, nan),
                 ('5', 'eg', 'EG', '1', nan, 12, nan, 0.3333333333333333, 179104, nan),
                 ('1', ROOT, 'Root', '', 6, 36, 1.0, 1.0, 1000000, nan),
                 ('37', 'g', 'G', '35', nan, 7, nan, 0.19444444444444445, 78575, nan),
                 ('33', 'd', 'D', '31', nan, 4, nan, 0.1111111111111111, 44900, nan),
                 ('31', 'cde', '... CDE', '28', 3, 12, 0.5, 0.3333333333333333, 134701, nan),
                 ('19', 'cf', 'CF', '11', 3, 9, 0.5, 0.25, 115138, nan),
                 ('28', 'cde+++', 'cde+++ ...', '27', 3, 12, 0.5, 0.3333333333333333, 134701, nan),
                 ('35', 'eg', 'EG', '27', nan, 12, nan, 0.3333333333333333, 134701, nan),
                 ('17', 'd', 'D', '15', nan, 4, nan, 0.1111111111111111, 51172, nan),
                 ('2', 'ab', 'AB', '1', 3, 3, 0.5, 0.08333333333333333, 44776, nan),
                 ('24', 'cdeeg+++', 'cdeeg+++ ...', '22', 3, 19, 0.5, 0.5277777777777778, 269402,
                  nan),
                 ('4', 'b', 'B', '2', 2, 2, 0.3333333333333333, 0.05555555555555555, 29850, nan)}
        lines = set(data.get_col())
        # from ontosunburst.sunburst_fig import generate_sunburst_fig
        # generate_sunburst_fig(data, 'test')
        self.assertEqual(len(lines), len(exp_l))
        for line in lines:
            line = tuple([nan if type(x) != str and np.isnan(x) else x for x in line])
            self.assertIn(line, exp_l)
