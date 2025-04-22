from ontosunburst.ontosunburst import *
import argparse


def get_command_line_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, required=True, help='Interest set input')
    parser.add_argument('--ref', '-r', type=str, required=False, help='Reference set input')
    parser.add_argument('--ontology', '--onto', type=str, required=False, help='Ontology used')
    parser.add_argument('--input_root', '-ir', type=str, required=False, help='Ontology root')
    parser.add_argument('--analysis', '-a', type=str, required=False, default=TOPOLOGY_A,
                        help='Type of analysis')
    parser.add_argument('--output', '-o', type=str, required=False, default='sunburst',
                        help='Output path+name')
    parser.add_argument('--ontology_dag', '-od', type=str, required=False,
                        help='Class ontology json file')
    parser.add_argument('--id_to_labels', '-itl', type=str, required=False,
                        help='IDs to labels json file')
    parser.add_argument('--no_labels', '-nl', action='store_false', required=False, default=True,
                        help='Show labels ')
    parser.add_argument('--test', '-t', type=str, required=False, default=BINOMIAL_TEST,
                        help='Enrichment stat test')
    parser.add_argument('--rcut', type=str, required=False, default=ROOT_CUT,
                        help='Type of root cut')
    parser.add_argument('--pcut', type=str, required=False, default=PATH_UNCUT,
                        help='Type of path cut')
    parser.add_argument('--r_base', action='store_true', required=False, default=False,
                        help='Reference base')
    parser.add_argument('--show_leaves', '-sl',  action='store_true', default=False, required=False,
                        help='Show leaves')
    parser.add_argument('--kwargs', nargs=argparse.REMAINDER, help="Additional keyword arguments")
    args = parser.parse_args()
    return args


def main():
    args = get_command_line_args()
    kwargs = {}
    if args.kwargs:
        for arg in args.kwargs:
            if "=" in arg:
                key, value = arg.split("=", 1)
                kwargs[key] = value
            else:
                raise ValueError(f"Argument {arg} is not in the form key=value")
    reference_set, ref_abundances, scores = extract_input(args.ref)
    metabolic_objects, abundances, scores = extract_input(args.input)
    ontosunburst(interest_set=metabolic_objects,
                 ontology=args.ontology,
                 input_root=args.input_root,
                 abundances=abundances,
                 scores=scores,
                 reference_set=reference_set,
                 ref_abundances=ref_abundances,
                 analysis=args.analysis,
                 output=args.output,
                 write_output=True,
                 ontology_dag_input=args.ontology_dag,
                 id_to_label_input=args.id_to_labels,
                 labels=args.no_labels,
                 test=args.test,
                 root_cut=args.rcut,
                 path_cut=args.pcut,
                 ref_base=args.r_base,
                 show_leaves=args.show_leaves,
                 **kwargs)


def extract_input(input_file):
    if input_file is not None:
        id_lst = []
        ab_lst = []
        sc_d = {}
        with open(input_file, 'r') as f:
            for l in f:
                l = l.strip().split('\t')
                id_lst.append(l[0])
                if len(l) == 3:
                    sc_d[l[0]] = float(l[2])
                if len(l) >= 2:
                    ab_lst.append(float(l[1]))
                else:
                    ab_lst.append(1)
        if not sc_d:
            sc_d = None
        return id_lst, ab_lst, sc_d
    return None, None, None
