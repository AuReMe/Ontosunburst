import unittest
from ontosunburst.ontosunburst import *

# METACYC
# ==================================================================================================

MET_SET = {'CPD-24674', 'CPD-24687', 'CPD-24688'}
REF_MET = {'CPD-24674', 'CPD-24687', 'CPD-24688',
           'CPD-12782', 'CPD-12784', 'CPD-12787',
           'CPD-12788', 'CPD-12789', 'CPD-12796',
           'CPD-12797', 'CPD-12798', 'CPD-12805',
           'CPD-12806', 'CPD-12812', 'CPD-12816',
           'CPD-1282', 'CPD-12824', 'CPD-1283'}

RXN_SET = {'CROTCOALIG-RXN', 'CYSTHIOCYS-RXN', 'NQOR-RXN'}
REF_RXN = {'CROTCOALIG-RXN', 'CYSTHIOCYS-RXN', 'NQOR-RXN',
           'RXN-14859', 'RXN-14873', 'RXN-14920',
           'RXN-14939', 'RXN-14975', 'RXN-21632',
           'RXN-21638', 'RXN-21652', 'RXN-8954'}

PWY_SET = {'2ASDEG-PWY', '4AMINOBUTMETAB-PWY', 'ALLANTOINDEG-PWY'}
REF_PWY = {'2ASDEG-PWY', '4AMINOBUTMETAB-PWY', 'ALLANTOINDEG-PWY',
           'CRNFORCAT-PWY', 'PWY-7195', 'PWY-7219',
           'PWY-7251', 'PWY-7351', 'PWY-7401',
           'PWY18C3-22', 'PWY0-1600', 'SERDEG-PWY'}


class MetacycTest(unittest.TestCase):
    def test_cpd_metacyc_proportion(self):
        metacyc_ontosunburst(metabolic_objects=MET_SET, output='test')

    def test_cpd_metacyc_comparison(self):
        metacyc_ontosunburst(metabolic_objects=MET_SET, reference_set=REF_MET, output='test')

    def test_rxn_metacyc_proportion(self):
        metacyc_ontosunburst(metabolic_objects=RXN_SET, output='test')

    def test_rxn_metacyc_comparison(self):
        metacyc_ontosunburst(metabolic_objects=RXN_SET, reference_set=REF_RXN, output='test')

    def test_pwy_metacyc_proportion(self):
        metacyc_ontosunburst(metabolic_objects=PWY_SET, output='test')

    def test_pwy_metacyc_comparison(self):
        metacyc_ontosunburst(metabolic_objects=PWY_SET, reference_set=REF_PWY, output='test')


# EC
# ==================================================================================================

EC_SET = {'2.6.1.45', '1.1.1.25', '1.1.1.140'}
REF_EC = {'2.6.1.45', '1.1.1.25', '1.1.1.140',
          '1.14.14.52', '2.7.1.137', '7.1.1.8',
          '1.17.4.5', '2.3.1.165', '3.2.1.53',
          '3.2.1.91', '6.3.4.2', '5.4.99.8'}


class EcTest(unittest.TestCase):

    def test_ec_proportion(self):
        ec_ontosunburst(ec_set=REF_EC, output='test')

    def test_ec_comparison(self):
        ec_ontosunburst(ec_set=EC_SET, reference_set=REF_EC, output='test')
