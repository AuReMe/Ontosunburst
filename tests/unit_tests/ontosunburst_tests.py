import unittest
import io
from functools import wraps

from ontosunburst.ontosunburst import *
from ontosunburst.dag2tree import *

"""
Tests manually good file creation.
No automatic tests integrated.
"""

# ==================================================================================================
# GLOBAL
# ==================================================================================================

# Complex structure (multi-parental)
ROOT = 'root'
C_LST = ['a', 'b', 'c']
C_REF = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
C_LAB = [1, 2, 3]
C_RAB = [1, 2, 3, 4, 5, 6, 7, 8]
C_ONTO = {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf'], 'd': ['cde'], 'e': ['cde', 'eg'],
          'f': ['cf'], 'g': ['gh', 'eg'], 'h': ['gh'],
          'ab': [ROOT], 'cde': ['cdecf', 'cdeeg'], 'cf': ['cdecf'],
          'eg': [ROOT, 'cdeeg'], 'gh': [ROOT],
          'cdecf': [ROOT], 'cdeeg': ['cdeeg+'], 'cdeeg+': [ROOT]}
C_LABELS = {ROOT: 'Root', 'cdeeg+': 'CDEEG+', 'cdeeg': 'CDEEG', 'cdecf': 'CDECF', 'gh': 'GH',
            'eg': 'EG', 'cde': 'CDE', 'cf': 'CF', 'h': 'H', 'g': 'G', 'f': 'F', 'e': 'E', 'd': 'D',
            'c': 'C', 'ab': 'AB', 'b': 'B'}

# Enrichment values
E_LST = ['02', '03', '04', '05', '08', '09']
E_LAB = [23, 20, 1, 4, 1, 1]
E_REF = ['01', '02', '03', '04', '05', '06', '07', '08', '09']
E_RAB = [14, 26, 20, 10, 20, 5, 1, 1, 3]
E_ONTO = {'01': ['00'], '02': ['00'], '03': ['00'], '04': ['00'], '05': ['01'],
          '06': ['01'], '07': ['01'], '08': ['02'], '09': ['02']}
E_LABElS = {'00': '0', '01': '1', '02': '2', '03': '3', '04': '4',
            '05': '5', '06': '6', '07': '7', '08': '8', '09': '9'}

# GLOBAL VALUES
G_LST = ['a', 'b', 'c']
G_REF = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
G_LAB = [23, 20, 5]
G_RAB = [14, 26, 20, 10, 20, 5, 4, 3, 1]
G_ONTO = {'a': ['i'], 'b': ['i'], 'c': ['j', 'k'], 'd': ['j'], 'e': ['j', 'l'],
          'f': ['k'], 'g': ['m', 'l'], 'h': ['m'],
          'i': ['x'], 'j': ['n', 'o'], 'k': ['n'],
          'l': ['x', 'o'], 'm': ['x'],
          'n': ['x'], 'o': ['v'], 'v': ['w'], 'w': ['x'], 'x': ['r'],
          'p': ['i'], 'q': ['p'], 's': ['q'], 't': ['p'], 'u': ['t']}
G_LABELS = {'r': 'Root', 'v': 'V', 'o': 'O', 'n': 'N', 'm': 'M',
            'l': 'L', 'j': 'J', 'k': 'K', 'h': 'H', 'g': 'G', 'f': 'F', 'e': 'E', 'd': 'D',
            'c': 'C', 'i': 'I', 'b': 'B'}


# ==================================================================================================
# FUNCTIONS UTILS
# ==================================================================================================

def save_fig_json(fig, file):
    fig = fig.to_dict()
    with open(file, 'w') as f:
        json.dump(fig, f)


def sort_fig_dict_lst(fig):
    """
    {data: [{hovertext: [str],
             ids: [str],
             labels: [str],
             marker: {colors: [int]},
             parents: [str],
             values: [int]}]}
    """
    lst_keys = ['hovertext', 'ids', 'labels', 'parents', 'values']
    for k in lst_keys:
        fig['data'][0][k] = sorted(fig['data'][0][k])
    fig['data'][0]['marker']['colors'] = sorted(fig['data'][0]['marker']['colors'])
    return fig


def fig_to_lines(fig_dict):
    lst_keys = ['hovertext', 'ids', 'labels', 'parents', 'values']
    lines = set()
    for i in range(len(fig_dict['data'][0]['ids'])):
        line = tuple(fig_dict['data'][0][k][i] for k in lst_keys)
        line = line + (str(fig_dict['data'][0]['marker']['colors'][i]),)
        lines.add(line)
    return lines


def are_fig_dict_equals(fig1, fig2_file):
    fig1 = fig1.to_dict()
    fig1_l = fig_to_lines(fig1)
    fig1 = json.dumps(sort_fig_dict_lst(fig1), sort_keys=True)
    with open(fig2_file, 'r') as f:
        fig2 = json.load(f)
        fig2_l = fig_to_lines(fig2)
        fig2 = json.dumps(sort_fig_dict_lst(fig2), sort_keys=True)
    return (fig1 == fig2) and (fig1_l == fig2_l)


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

# Write image : fig.write_image('test_go_prop.png', width=1900, height=1000, scale=1)


# ==================================================================================================
# UNIT TESTS
# ==================================================================================================

# CUSTOM ONTO
# --------------------------------------------------------------------------------------------------

class TestOntosunburstCustomOnto(unittest.TestCase):

    def test_art(self):
        ontosunburst(interest_set=G_LST, input_root='r', ontology_dag_input=G_ONTO,
                     reference_set=None, abundances=None, ref_abundances=None, labels=True,
                     write_output=False, output='test_art', show_leaves=True, root_cut='uncut',
                     bg_color='white', ref_base=False, analysis=TOPOLOGY_A, path_cut=PATH_UNCUT)

    # TOPOLOGY : CUSTOM ONTO
    @test_for(ontosunburst)
    def test_ontosunburst_1(self):
        fig = ontosunburst(interest_set=C_LST, ontology=None, input_root=ROOT,
                           abundances=C_LAB, reference_set=C_REF, ref_abundances=C_RAB,
                           analysis='topology', output='test1', write_output=False,
                           ontology_dag_input=C_ONTO, id_to_label_input=C_LABELS, labels=True,
                           test=BINOMIAL_TEST, root_cut=ROOT_TOTAL_CUT, ref_base=True,
                           show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test1.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

    @test_for(ontosunburst)
    def test_ontosunburst_2(self):
        fig = ontosunburst(interest_set=C_REF, ontology=None, input_root=ROOT,
                           abundances=C_RAB, reference_set=None, ref_abundances=None,
                           analysis='topology', output='test2', write_output=False,
                           ontology_dag_input=C_ONTO, id_to_label_input=C_LABELS, labels=False,
                           test=BINOMIAL_TEST, root_cut=ROOT_CUT, ref_base=False, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test2.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

    @test_for(ontosunburst)
    def test_ontosunburst_3(self):
        fig = ontosunburst(interest_set=C_LST, ontology=None, input_root=ROOT,
                           abundances=C_LAB, reference_set=C_REF, ref_abundances=C_RAB,
                           analysis='topology', output='test3', write_output=False,
                           ontology_dag_input=C_ONTO, id_to_label_input=C_LABELS, labels=True,
                           test=BINOMIAL_TEST, root_cut=ROOT_UNCUT, ref_base=False,
                           show_leaves=False, bg_color='#ddddcc')
        w_fig_file = os.path.join('test_files', 'test3.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

    # ENRICHMENT : CUSTOM ONTO

    @test_for(ontosunburst)
    def test_ontosunburst_4(self):
        fig = ontosunburst(interest_set=E_LST, ontology=None, input_root='00',
                           abundances=E_LAB, reference_set=E_REF, ref_abundances=E_RAB,
                           analysis='enrichment', output='test4', write_output=False,
                           ontology_dag_input=E_ONTO, id_to_label_input=E_LABElS, labels=True,
                           test=BINOMIAL_TEST, root_cut=ROOT_CUT, ref_base=True, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test4.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

    @test_for(ontosunburst)
    def test_ontosunburst_5(self):
        fig = ontosunburst(interest_set=E_LST, ontology=None, input_root='00',
                           abundances=E_LAB, reference_set=E_REF, ref_abundances=E_RAB,
                           analysis='enrichment', output='test5', write_output=False,
                           ontology_dag_input=E_ONTO, id_to_label_input=E_LABElS, labels=True,
                           test=HYPERGEO_TEST, root_cut=ROOT_UNCUT, ref_base=False,
                           show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test5.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

    @test_for(ontosunburst)
    def test_ontosunburst_6(self):
        scores = {'01': 0.002, '02': 0.00007, '03': 0.9, '05': 0.004, '08': 0.2, '09': 0.0000012}
        fig = ontosunburst(interest_set=E_LST, ontology=None, input_root='00', scores=scores,
                           abundances=E_LAB, reference_set=E_REF, ref_abundances=E_RAB,
                           analysis='enrichment', output='test6', write_output=False,
                           ontology_dag_input=E_ONTO, id_to_label_input=E_LABElS, labels=True,
                           root_cut=ROOT_UNCUT, ref_base=False, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test6.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))


# METACYC
# --------------------------------------------------------------------------------------------------

MET_LST = ['CPD-24674', 'CPD-24687', 'CPD-24688']
MET_REF = ['CPD-24674', 'CPD-24687', 'CPD-24688',
           'CPD-12782', 'CPD-12784', 'CPD-12787',
           'CPD-12788', 'CPD-12789', 'CPD-12796',
           'CPD-12797', 'CPD-12798', 'CPD-12805',
           'CPD-12806', 'CPD-12812', 'CPD-12816',
           'CPD-1282', 'CPD-12824', 'CPD-1283']
MET_LAB = [1, 2, 3]
MET_RAB = [1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9]

RXN_LST = ['CROTCOALIG-RXN', 'CYSTHIOCYS-RXN', 'NQOR-RXN']
RXN_REF = ['CROTCOALIG-RXN', 'CYSTHIOCYS-RXN', 'NQOR-RXN',
           'RXN-14859', 'RXN-14873', 'RXN-14920',
           'RXN-14939', 'RXN-14975', 'RXN-21632',
           'RXN-21638', 'RXN-21652', 'RXN-8954']

PWY_LST = ['2ASDEG-PWY', '4AMINOBUTMETAB-PWY', 'ALLANTOINDEG-PWY']
PWY_REF = ['2ASDEG-PWY', '4AMINOBUTMETAB-PWY', 'ALLANTOINDEG-PWY',
           'CRNFORCAT-PWY', 'PWY-7195', 'PWY-7219',
           'PWY-7251', 'PWY-7351', 'PWY-7401',
           'PWY18C3-22', 'PWY0-1600', 'SERDEG-PWY']


class TestOntosunburstMetaCyc(unittest.TestCase):

    # Topology
    @test_for(ontosunburst)
    def test_ontosunburst_mc1(self):
        fig = ontosunburst(interest_set=MET_LST, ontology=METACYC,
                           abundances=MET_LAB, reference_set=MET_REF, ref_abundances=MET_RAB,
                           analysis=TOPOLOGY_A, output='mc1', write_output=False,
                           ontology_dag_input=None, labels=True,
                           root_cut=ROOT_CUT, ref_base=True, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test_mc1.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

    @test_for(ontosunburst)
    def test_ontosunburst_mc2(self):
        fig = ontosunburst(interest_set=MET_LST, ontology=METACYC,
                           abundances=MET_LAB, reference_set=MET_REF, ref_abundances=MET_RAB,
                           analysis=TOPOLOGY_A, output='mc2', write_output=False,
                           ontology_dag_input=None, id_to_label_input={'Acids': 'ACIDS !!'},
                           test=HYPERGEO_TEST, root_cut=ROOT_CUT, ref_base=False, show_leaves=False)
        w_fig_file = os.path.join('test_files', 'test_mc2.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

    @test_for(ontosunburst)
    def test_ontosunburst_mc3(self):
        fig = ontosunburst(interest_set=MET_REF, ontology=METACYC,
                           abundances=None, reference_set=None, ref_abundances=None,
                           analysis=TOPOLOGY_A, output='mc3', write_output=False,
                           ontology_dag_input=None, test=HYPERGEO_TEST, root_cut=ROOT_CUT,
                           ref_base=False, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test_mc3.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

    # Enrichment
    @test_for(ontosunburst)
    def test_ontosunburst_mc4(self):
        fig = ontosunburst(interest_set=MET_LST, ontology=METACYC,
                           abundances=MET_LAB, reference_set=MET_REF, ref_abundances=MET_RAB,
                           analysis=ENRICHMENT_A, output='mc4', write_output=False,
                           test=HYPERGEO_TEST, root_cut=ROOT_CUT, ref_base=False, show_leaves=False)
        w_fig_file = os.path.join('test_files', 'test_mc4.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

    @test_for(ontosunburst)
    def test_ontosunburst_mc5(self):
        fig = ontosunburst(interest_set=MET_LST, ontology=METACYC,
                           abundances=MET_LAB, reference_set=MET_REF, ref_abundances=MET_RAB,
                           analysis=ENRICHMENT_A, output='mc5', write_output=False,
                           test=BINOMIAL_TEST, root_cut=ROOT_CUT, ref_base=True, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test_mc5.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

    # Pathways
    @test_for(ontosunburst)
    def test_ontosunburst_mc6(self):
        fig = ontosunburst(interest_set=PWY_LST, ontology=METACYC,
                           abundances=None, reference_set=PWY_REF, ref_abundances=None,
                           analysis=TOPOLOGY_A, output='mc6', write_output=False,
                           test=BINOMIAL_TEST, root_cut=ROOT_CUT,
                           ref_base=True, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test_mc6.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

    # Reactions
    @test_for(ontosunburst)
    def test_ontosunburst_mc7(self):
        fig = ontosunburst(interest_set=RXN_LST, ontology=METACYC,
                           abundances=None, reference_set=RXN_REF, ref_abundances=None,
                           analysis=TOPOLOGY_A, output='mc7', write_output=False,
                           test=BINOMIAL_TEST, root_cut=ROOT_CUT,
                           ref_base=True, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test_mc7.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

    @test_for(ontosunburst)
    def test_ontosunburst_mc8(self):
        fig = ontosunburst(interest_set=MET_LST, ontology=METACYC,
                           abundances=MET_LAB, reference_set=MET_REF, ref_abundances=MET_RAB,
                           analysis=TOPOLOGY_A, output='mc8', write_output=False,
                           test=HYPERGEO_TEST, root_cut=ROOT_CUT, path_cut=PATH_BOUND,
                           ref_base=False, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test_mc8.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

    @test_for(ontosunburst)
    def test_ontosunburst_mc9(self):
        fig = ontosunburst(interest_set=MET_LST, ontology=METACYC,
                           abundances=MET_LAB, reference_set=MET_REF, ref_abundances=MET_RAB,
                           analysis=TOPOLOGY_A, output='mc9', write_output=False,
                           test=HYPERGEO_TEST, root_cut=ROOT_CUT, path_cut=PATH_HIGHER,
                           ref_base=True, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test_mc9.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))


# EC
# --------------------------------------------------------------------------------------------------

EC_LST = ['2.6.1.45', '1.1.1.25', '1.1.1.140', '1.1.2.-']
REF_EC = ['2.6.1.45', '1.1.1.25', '1.1.1.140',
          '1.14.14.52', '2.7.1.137', '7.1.1.8',
          '1.17.4.5', '2.3.1.165', '3.2.1.53',
          '3.2.1.91', '6.3.4.2', '5.4.99.8']


class TestOntosunburstEC(unittest.TestCase):

    @test_for(ontosunburst)
    def test_ontosunburst_ec1(self):
        fig = ontosunburst(interest_set=EC_LST, ontology=EC,
                           abundances=None, reference_set=REF_EC, ref_abundances=None,
                           analysis=TOPOLOGY_A, output='ec1', write_output=False,
                           test=BINOMIAL_TEST, root_cut=ROOT_CUT,
                           ref_base=True, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test_ec1.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

    @test_for(ontosunburst)
    def test_ontosunburst_ec2(self):
        fig = ontosunburst(interest_set=EC_LST, ontology=EC,
                           abundances=None, reference_set=REF_EC, ref_abundances=None,
                           analysis=TOPOLOGY_A, output='ec2', write_output=False,
                           labels=False, test=BINOMIAL_TEST, root_cut=ROOT_CUT,
                           ref_base=True, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test_ec2.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

    @test_for(ontosunburst)
    def test_ontosunburst_ec3(self):
        fig = ontosunburst(interest_set=EC_LST, ontology=EC,
                           abundances=None, reference_set=REF_EC, ref_abundances=None,
                           analysis=ENRICHMENT_A, output='ec3', write_output=False,
                           class_ontology=None, id_to_label_input={'1.-.-.-': ':D'},
                           test=BINOMIAL_TEST, root_cut=ROOT_CUT,
                           ref_base=True, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test_ec3.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))


# KEGG
# --------------------------------------------------------------------------------------------------

KEGG_EX = ['M00572', 'M00308', 'M00844', 'M00633', 'M00176', 'M00535', 'M00573', 'M00970', 'M00131',
           'M00620']


class TestOntosunburstKegg(unittest.TestCase):

    @test_for(ontosunburst)
    def test_ontosunburst_kg1(self):
        fig = ontosunburst(interest_set=KEGG_EX, ontology=KEGG,
                           abundances=None, reference_set=None, ref_abundances=None,
                           analysis=TOPOLOGY_A, output='kg1', write_output=False,
                           test=BINOMIAL_TEST, root_cut=ROOT_CUT,
                           ref_base=False, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test_kg1.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))


# CHEBI / CHEBI ROLES
# --------------------------------------------------------------------------------------------------

CH_LST = ['CHEBI:38028', 'CHEBI:28604', 'CHEBI:85146']
REF_CH = ['CHEBI:38028', 'CHEBI:28604', 'CHEBI:85146',
          'CHEBI:23066', 'CHEBI:27803', 'CHEBI:37565',
          'CHEBI:58215', 'CHEBI:79983', 'CHEBI:42639']


class TestOntosunburstChEBI(unittest.TestCase):

    @test_for(ontosunburst)
    def test_ontosunburst_ch1(self):
        fig = ontosunburst(interest_set=CH_LST, ontology=CHEBI, reference_set=REF_CH,
                           abundances=None, ref_abundances=None,
                           analysis=TOPOLOGY_A, output='ch1', write_output=False,
                           test=BINOMIAL_TEST, root_cut=ROOT_CUT,
                           ref_base=True, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test_ch1.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

    @test_for(ontosunburst)
    def test_ontosunburst_ch2(self):
        fig = ontosunburst(interest_set=CH_LST, ontology=CHEBI, reference_set=REF_CH,
                           abundances=None, ref_abundances=None, labels=False,
                           analysis=TOPOLOGY_A, output='ch2', write_output=False,
                           test=BINOMIAL_TEST, root_cut=ROOT_CUT,
                           ref_base=False, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test_ch2.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

    @test_for(ontosunburst)
    def test_ontosunburst_ch3(self):
        fig = ontosunburst(interest_set=CH_LST, ontology=CHEBI, reference_set=REF_CH,
                           abundances=None, ref_abundances=None, labels=True,
                           analysis=ENRICHMENT_A, output='ch3', write_output=False,
                           test=BINOMIAL_TEST, root_cut=ROOT_CUT,
                           ref_base=True, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test_ch3.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))


class TestOntosunburstChEBIRoles(unittest.TestCase):

    @test_for(ontosunburst)
    def test_ontosunburst_chr1(self):
        fig = ontosunburst(interest_set=CH_LST, ontology=CHEBI_R, reference_set=REF_CH,
                           abundances=None, ref_abundances=None,
                           analysis=TOPOLOGY_A, output='chr1', write_output=False,
                           test=BINOMIAL_TEST, root_cut=ROOT_CUT,
                           ref_base=True, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test_chr1.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

    @test_for(ontosunburst)
    def test_ontosunburst_chr2(self):
        fig = ontosunburst(interest_set=CH_LST, ontology=CHEBI_R, reference_set=REF_CH,
                           abundances=None, ref_abundances=None, labels=False,
                           analysis=TOPOLOGY_A, output='chr2', write_output=False,
                           test=BINOMIAL_TEST, root_cut=ROOT_CUT,
                           ref_base=False, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test_chr2.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

    @test_for(ontosunburst)
    def test_ontosunburst_chr3(self):
        fig = ontosunburst(interest_set=CH_LST, ontology=CHEBI_R, reference_set=REF_CH,
                           abundances=None, ref_abundances=None, labels=True,
                           analysis=ENRICHMENT_A, output='chr3', write_output=False,
                           test=BINOMIAL_TEST, root_cut=ROOT_CUT,
                           ref_base=True, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test_chr3.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))


# GO / GO_BP / GO_MF / GO_CC
# --------------------------------------------------------------------------------------------------

GO_BP_LST = ['GO:0061948', 'GO:1902266', 'GO:0003046']
GO_MF_LST = ['GO:0051870', 'GO:0010013', 'GO:0031386']
GO_CC_LST = ['GO:0098061', 'GO:1990177', 'GO:0098833']
GO_LST = ['GO:0061948', 'GO:1902266', 'GO:0003046',
          'GO:0051870', 'GO:0010013', 'GO:0031386',
          'GO:0098061', 'GO:1990177', 'GO:0098833']


class TestOntosunburstGO(unittest.TestCase):

    @test_for(ontosunburst)
    def test_ontosunburst_go_bp1(self):
        fig = ontosunburst(interest_set=GO_BP_LST, ontology=GO_BP, reference_set=None,
                           abundances=None, ref_abundances=None,
                           analysis=TOPOLOGY_A, output='go_bp1', write_output=False,
                           test=BINOMIAL_TEST, root_cut=ROOT_CUT,
                           ref_base=False, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test_go_bp1.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

    @test_for(ontosunburst)
    def test_ontosunburst_go_mf1(self):
        fig = ontosunburst(interest_set=GO_MF_LST, ontology=GO_MF, reference_set=None,
                           abundances=None, ref_abundances=None,
                           analysis=TOPOLOGY_A, output='go_mf1', write_output=False,
                           test=BINOMIAL_TEST, root_cut=ROOT_CUT,
                           ref_base=False, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test_go_mf1.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

    @test_for(ontosunburst)
    def test_ontosunburst_go_cc1(self):
        fig = ontosunburst(interest_set=GO_CC_LST, ontology=GO_CC, reference_set=None,
                           abundances=None, ref_abundances=None,
                           analysis=TOPOLOGY_A, output='go_cc1', write_output=False,
                           test=BINOMIAL_TEST, root_cut=ROOT_CUT,
                           ref_base=False, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test_go_cc1.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

    @test_for(ontosunburst)
    def test_ontosunburst_go1(self):
        fig = ontosunburst(interest_set=GO_LST, ontology=GO, reference_set=None,
                           abundances=None, ref_abundances=None,
                           analysis=TOPOLOGY_A, output='go1', write_output=False,
                           test=BINOMIAL_TEST, root_cut=ROOT_CUT,
                           ref_base=False, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test_go1.json')
        # save_fig_json(fig, w_fig_file)
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

