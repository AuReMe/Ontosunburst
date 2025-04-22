import json
import os.path
import unittest
import io
import sys
import copy
from functools import wraps
from unittest.mock import patch

from ontosunburst.tree2sunburst import *

"""
Tests manually good file creation.
No automatic tests integrated.
"""

# ==================================================================================================
# GLOBAL
# ==================================================================================================

# Complex topology
# --------------------------------------------------------------------------------------------------

ROOT = 'root'
MC_ONTO = {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf'], 'd': ['cde'], 'e': ['cde', 'eg'],
           'f': ['cf'], 'g': ['gh', 'eg'], 'h': ['gh'],
           'ab': [ROOT], 'cde': ['cdecf', 'cdeeg'], 'cf': ['cdecf'],
           'eg': ['cdeeg', ROOT], 'gh': [ROOT],
           'cdecf': [ROOT], 'cdeeg': ['cdeeg+'], 'cdeeg+': [ROOT]}

MC_AB = {ROOT: 6, 'cde': 3, 'cf': 3, 'cdecf': 3, 'cdeeg+': 3, 'cdeeg': 3, 'c': 3, 'ab': 3,
         'b': 2, 'a': 1}
MC_REF_AB = {ROOT: 36, 'cdeeg+': 19, 'cdeeg': 19, 'cdecf': 18, 'gh': 15, 'eg': 12, 'cde': 12,
             'cf': 9, 'h': 8, 'g': 7, 'f': 6, 'e': 5, 'd': 4, 'c': 3, 'ab': 3, 'b': 2, 'a': 1}

MC_LABELS = {ROOT: 'Root', 'cdeeg+': 'CDEEG+', 'cdeeg': 'CDEEG', 'cdecf': 'CDECF', 'gh': 'GH',
             'eg': 'EG', 'cde': 'CDE', 'cf': 'CF', 'h': 'H', 'g': 'G', 'f': 'F', 'e': 'E', 'd': 'D',
             'c': 'C', 'ab': 'AB', 'b': 'B'}
MC_DATA = TreeData()
MC_DATA.dag_to_tree(MC_AB, MC_REF_AB, MC_ONTO, ROOT, MC_LABELS)
MC_DATA.calculate_proportions(True)

# Enrichment
E_AB = {'00': 50, '01': 5, '02': 25, '03': 20, '04': 1, '05': 5, '06': nan, '07': nan,
        '08': 1, '09': 1}
E_REF_AB = {'00': 100, '01': 40, '02': 30, '03': 20, '04': 10, '05': 20, '06': 5, '07': 1,
            '08': 1, '09': 3}
E_LABElS = {'00': '0', '01': '1', '02': '2', '03': '3', '04': '4',
            '05': '5', '06': '6', '07': '7', '08': '8', '09': '9'}
E_ONTO = {'01': ['00'], '02': ['00'], '03': ['00'], '04': ['00'], '05': ['01'],
          '06': ['01'], '07': ['01'], '08': ['02'], '09': ['02']}

E_DATA = TreeData()
E_DATA.dag_to_tree(E_AB, E_REF_AB, E_ONTO, '00', E_LABElS, True)
E_DATA.calculate_proportions(True)
E_SIGN = E_DATA.make_enrichment_analysis(BINOMIAL_TEST)


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


# ==================================================================================================
# UNIT TESTS
# ==================================================================================================

# TEST
# --------------------------------------------------------------------------------------------------

class TestSunburstFigure(unittest.TestCase):
    @test_for(get_fig_kwargs)
    def test_get_kwargs_topology(self):
        c_min, c_max, c_mid, max_depth, colorscale, title, colorbar_legend, background_color, \
            font_color, font_size, table_title, table_legend, table_color = \
            get_fig_kwargs('out', TOPOLOGY_A)
        self.assertEqual(c_min, 1)
        self.assertEqual(c_max, None)
        self.assertEqual(c_mid, None)
        self.assertEqual(max_depth, 7)
        self.assertEqual(colorscale, px.colors.get_colorscale('Viridis'))
        self.assertEqual(title, 'out : Proportion of classes')
        self.assertEqual(colorbar_legend, 'Count')
        self.assertEqual(background_color, 'rgba(255, 255, 255, 0)')
        self.assertEqual(font_color, '#111111')
        self.assertEqual(font_size, 20)
        self.assertEqual(table_title, 'Significant p-values')
        self.assertEqual(table_legend, 'IDs')
        self.assertEqual(table_color, '#666666')

    @test_for(get_fig_kwargs)
    def test_get_kwargs_enrichment(self):
        c_min, c_max, c_mid, max_depth, colorscale, title, colorbar_legend, background_color, \
            font_color, font_size, table_title, table_legend, table_color = \
            get_fig_kwargs('out', ENRICHMENT_A)
        self.assertEqual(c_min, -10)
        self.assertEqual(c_max, 10)
        self.assertEqual(c_mid, 0)
        self.assertEqual(max_depth, 7)
        self.assertEqual(colorscale, px.colors.get_colorscale('RdBu'))
        self.assertEqual(title, 'out : Classes enrichment representation')
        self.assertEqual(colorbar_legend, 'Log10(p-value)')
        self.assertEqual(background_color, 'rgba(255, 255, 255, 0)')
        self.assertEqual(font_color, '#111111')
        self.assertEqual(font_size, 20)
        self.assertEqual(table_title, 'Significant p-values')
        self.assertEqual(table_legend, 'IDs')
        self.assertEqual(table_color, '#666666')

    @test_for(get_fig_kwargs)
    def test_get_kwargs(self):
        c_min, c_max, c_mid, max_depth, colorscale, title, colorbar_legend, background_color, \
            font_color, font_size, table_title, table_legend, table_color = \
            get_fig_kwargs('out', TOPOLOGY_A, c_min=4, c_max=8, c_mid=6, max_depth=10,
                           colorscale='Twilight', title='My title', colorbar_legend='Total',
                           bg_color='#000000', font_color='#ffffff', font_size=24,
                           table_title='Significant classes', table_legend='Name',
                           table_color='#222222')
        self.assertEqual(c_min, 4)
        self.assertEqual(c_max, 8)
        self.assertEqual(c_mid, 6)
        self.assertEqual(max_depth, 10)
        self.assertEqual(colorscale, px.colors.get_colorscale('Twilight'))
        self.assertEqual(title, 'My title')
        self.assertEqual(colorbar_legend, 'Total')
        self.assertEqual(background_color, '#000000')
        self.assertEqual(font_color, '#ffffff')
        self.assertEqual(font_size, 24)
        self.assertEqual(table_title, 'Significant classes')
        self.assertEqual(table_legend, 'Name')
        self.assertEqual(table_color, '#222222')

    @test_for(check_kwargs)
    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_check_kwargs(self, mock_stdout):
        check_kwargs(bg_color='black', cmin=2, max_detph=3, c_mif=10, framboise='fruit',
                     c_max='blue', backgroung_color='white')
        output = mock_stdout.getvalue().strip()
        expected_msg = 'Unknown kwarg "cmin", did you mean "c_min" ?\n' \
                       'Unknown kwarg "max_detph", did you mean "max_depth" ?\n' \
                       'Unknown kwarg "c_mif", did you mean "c_min" ?\n' \
                       'Unknown kwarg "framboise"\n' \
                       '"c_max" must be of type "<class \'float\'>" not "<class \'str\'>"\n' \
                       'Unknown kwarg "backgroung_color", did you mean "bg_color" ?'
        self.assertEqual(output, expected_msg)

    @test_for(get_hover_fig_text)
    def test_get_hover_fig_text_enrich_ref(self):
        data = copy.deepcopy(E_DATA)
        text_list = get_hover_fig_text(data, ENRICHMENT_A, True)
        self.assertEqual(len(text_list), 10)
        self.assertEqual(text_list[0], 'P value: 1.0<br>Weight: <b>50</b>'
                                       '<br>Reference weight: 100<br>Proportion: <b>100.0%</b>'
                                       '<br>Reference proportion: 100.0%<br>ID: 00')
        self.assertEqual(text_list[2], 'P value: 0.07883064215278136<br>Weight: <b>5</b>'
                                       '<br>Reference weight: 20<br>Proportion: <b>10.0%</b>'
                                       '<br>Reference proportion: 20.0%<br>ID: 05')

    @test_for(get_hover_fig_text)
    def test_get_hover_fig_text_enrich_no_ref(self):
        data = copy.deepcopy(E_DATA)
        text_list = get_hover_fig_text(data, ENRICHMENT_A, False)
        self.assertEqual(len(text_list), 10)
        self.assertEqual(text_list[0], 'P value: 1.0<br>Weight: <b>50</b><br>'
                                       'Reference weight: 100<br>Proportion: <b>100.0%</b><br>'
                                       'Reference proportion: 100.0%<br>ID: 00')
        self.assertEqual(text_list[2], 'P value: 0.07883064215278136<br>Weight: <b>5</b><br>'
                                       'Reference weight: 20<br>Proportion: <b>10.0%</b><br>'
                                       'Reference proportion: 20.0%<br>ID: 05')

    @test_for(get_hover_fig_text)
    def test_get_hover_fig_text_topology_ref(self):
        data = copy.deepcopy(E_DATA)
        text_list = get_hover_fig_text(data, TOPOLOGY_A, True)
        self.assertEqual(len(text_list), 10)
        self.assertEqual(text_list[0], 'Weight: <b>50</b><br>Reference weight: 100'
                                       '<br>Proportion: <b>100.0%</b>'
                                       '<br>Reference proportion: 100.0%<br>ID: 00')
        self.assertEqual(text_list[2], 'Weight: <b>5</b><br>Reference weight: 20'
                                       '<br>Proportion: <b>10.0%</b>'
                                       '<br>Reference proportion: 20.0%<br>ID: 05')

    @test_for(get_hover_fig_text)
    def test_get_hover_fig_text_topology_no_ref(self):
        data = copy.deepcopy(E_DATA)
        text_list = get_hover_fig_text(data, TOPOLOGY_A, False)
        self.assertEqual(len(text_list), 10)
        self.assertEqual(text_list[0], 'Weight: <b>50</b><br>Proportion: <b>100.0%</b><br>ID: 00')
        self.assertEqual(text_list[2], 'Weight: <b>5</b><br>Proportion: <b>10.0%</b><br>ID: 05')

    @test_for(generate_sunburst_fig)
    def test_generate_sunburst_fig_case1(self):
        data = copy.deepcopy(E_DATA)
        fig = generate_sunburst_fig(data, 'case1', analysis=ENRICHMENT_A, write_fig=False,
                                    test=HYPERGEO_TEST, ref_set=True, significant=E_SIGN)
        w_fig_file = os.path.join('test_files', 'fig_case1.json')
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

    @test_for(generate_sunburst_fig)
    def test_generate_sunburst_fig_case2(self):
        data = copy.deepcopy(E_DATA)
        fig = generate_sunburst_fig(data, 'case2', analysis=ENRICHMENT_A, write_fig=False,
                                    test=BINOMIAL_TEST, significant=E_SIGN)
        w_fig_file = os.path.join('test_files', 'fig_case2.json')
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

    @test_for(generate_sunburst_fig)
    def test_generate_sunburst_fig_case3(self):
        data = copy.deepcopy(E_DATA)
        fig = generate_sunburst_fig(data, 'case3', analysis=ENRICHMENT_A,
                                    test=HYPERGEO_TEST, significant=E_SIGN,
                                    root_cut=ROOT_UNCUT, write_fig=False,
                                    title='Another title', colorscale='PuOr_r',
                                    bg_color='#222222', font_color='#eeeeee', font_size=25,
                                    table_legend='Number')
        w_fig_file = os.path.join('test_files', 'fig_case3.json')
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))

    @test_for(generate_sunburst_fig)
    def test_generate_sunburst_fig_case4(self):
        data = copy.deepcopy(MC_DATA)
        fig = generate_sunburst_fig(data, 'case4', analysis=TOPOLOGY_A, bg_color='black',
                                    font_color='white', write_fig=False)
        w_fig_file = os.path.join('test_files', 'fig_case4.json')
        self.assertTrue(are_fig_dict_equals(fig, w_fig_file))
