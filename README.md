[![](https://img.shields.io/badge/python-3.10-blue.svg)]()
[![](https://img.shields.io/badge/version-0.0.5-green.svg)](https://github.com/AuReMe/Ontosunburst/releases/tag/v0.0.5)
[![](https://img.shields.io/badge/documentation-Wiki-orange.svg)](https://github.com/AuReMe/Ontosunburst/wiki)


# Ontosunburst

Sunburst visualisation of an ontology representing classes of sets
of metabolic objects


![image](./Figures/main_fig_topo.png)
![image](./Figures/main_fig_enrich.png)

## Requirements

### Mandatory
Python 3.10 recommended

Requirements from `requirements.txt`

- numpy>=1.26.1
- plotly>=5.17.0
- scipy>=1.11.3
- SPARQLWrapper>=2.0.0
- pandas>=1.5.3

### Optional

Need *Apache Jena Fuseki* SPARQL server to generate your own CHEBI and GO ontology input files

- Download *Apache Jena Fuseki* : https://jena.apache.org/download/index.cgi 
- Download ChEBI ontology : https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/
  (chebi.owl or chebi_lite.owl)
- Download GO ontology : https://geneontology.org/docs/download-ontology/ (go-basic.owl)

## Installation

### PyPI

```commandline
pip install ontosunburst
```

### Local

Inside the cloned repository :

```commandline
pip install -e .
```



## Utilisation

### Availabilities

#### 5 **Ontologies :**

With local files :
- MetaCyc (compounds, reactions, pathways)
- EC (EC-numbers)
- KEGG Ontology (modules, pathways, ko, ko_transporter, metabolite, metabolite_lipid)
- ChEBI (chebi: chebi classes + chebi_r: chebi roles)
- Gene Ontology (go_bp: biological process + go_mf: molecular function + go_cc:
cellular component + go: aggregation of 3)

Personal ontology possible :
- Define all the ontology classes relationship in 
a dictionary `{class: [parent classes]}`
- Define the root : unique class with no parents

#### 2 **Analysis :**

- Topology (**1 set** + 1 optional reference set) : displays proportion 
(number of occurrences) representation of all classes
- Enrichment (**1 set** + **1 reference set**) :  displays enrichment 
analysis significance of a set according to a reference set of metabolic 
objects

# Documentation

View full documentation here : https://github.com/AuReMe/Ontosunburst/wiki 

# References:
 
- **plotly:** *Plotly Technologies Inc. Collaborative data science. Montréal, QC, 2015. https://plot.ly.*

### Ontologies DataBases

- **MetaCyc:** 
  - Caspi, R., Billington, R., Keseler, I. M., Kothari, A., Krummenacker, M., Midford, P. E., 
  Ong, W. K., Paley, S., Subhraveti, P., & Karp, P. D. (2020). The MetaCyc database of metabolic 
  pathways and enzymes - a 2019 update. Nucleic acids research, 48(D1), D445–D453. 
  https://doi.org/10.1093/nar/gkz862
  - Ron Caspi, Kate Dreher, Peter D. Karp, The challenge of constructing, classifying, and 
  representing metabolic pathways, FEMS Microbiology Letters, Volume 345, Issue 2, August 2013, 
  Pages 85–93, https://doi.org/10.1111/1574-6968.12194
  - Karp, P.D., Caspi, R. A survey of metabolic databases emphasizing the MetaCyc family. Arch 
  Toxicol 85, 1015–1033 (2011). https://doi.org/10.1007/s00204-011-0705-2

- **KEGG:**
  - Minoru Kanehisa, Miho Furumichi, Yoko Sato, Yuriko Matsuura, Mari Ishiguro-Watanabe, KEGG: 
  biological systems database as a model of the real world, Nucleic Acids Research, Volume 53, 
  Issue D1, 6 January 2025, Pages D672–D677, https://doi.org/10.1093/nar/gkae909
- **EC:**
  - Amos Bairoch, The ENZYME database in 2000, Nucleic Acids Research, Volume 28, Issue 1, 1 
  January 2000, Pages 304–305, https://doi.org/10.1093/nar/28.1.304
- **ChEBI:**
  - Janna Hastings, Gareth Owen, Adriano Dekker, Marcus Ennis, Namrata Kale, Venkatesh 
  Muthukrishnan, Steve Turner, Neil Swainston, Pedro Mendes, Christoph Steinbeck, ChEBI in 2016: 
  Improved services and an expanding collection of metabolites, Nucleic Acids Research, Volume 44, 
  Issue D1, 4 January 2016, Pages D1214–D1219, https://doi.org/10.1093/nar/gkv1031
- **GO:**
  - The Gene Ontology Consortium , Suzi A Aleksander, James Balhoff, Seth Carbon, J Michael Cherry, 
  Harold J Drabkin, Dustin Ebert, Marc Feuermann, Pascale Gaudet, Nomi L Harris, David P Hill, 
  Raymond Lee, Huaiyu Mi, Sierra Moxon, Christopher J Mungall, Anushya Muruganugan, Tremayne 
  Mushayahama, Paul W Sternberg, Paul D Thomas, Kimberly Van Auken, Jolene Ramsey, Deborah A 
  Siegele, Rex L Chisholm, Petra Fey, Maria Cristina Aspromonte, Maria Victoria Nugnes, Federica 
  Quaglia, Silvio Tosatto, Michelle Giglio, Suvarna Nadendla, Giulia Antonazzo, Helen Attrill, 
  Gil dos Santos, Steven Marygold, Victor Strelets, Christopher J Tabone, Jim Thurmond, Pinglei 
  Zhou, Saadullah H Ahmed, Praoparn Asanitthong, Diana Luna Buitrago, Meltem N Erdol, Matthew C 
  Gage, Mohamed Ali Kadhum, Kan Yan Chloe Li, Miao Long, Aleksandra Michalak, Angeline Pesala, 
  Armalya Pritazahra, Shirin C C Saverimuttu, Renzhi Su, Kate E Thurlow, Ruth C Lovering, Colin 
  Logie, Snezhana Oliferenko, Judith Blake, Karen Christie, Lori Corbani, Mary E Dolan, Harold J 
  Drabkin, David P Hill, Li Ni, Dmitry Sitnikov, Cynthia Smith, Alayne Cuzick, James Seager, 
  Laurel Cooper, Justin Elser, Pankaj Jaiswal, Parul Gupta, Pankaj Jaiswal, Sushma Naithani, 
  Manuel Lera-Ramirez, Kim Rutherford, Valerie Wood, Jeffrey L De Pons, Melinda R Dwinell, G 
  Thomas Hayman, Mary L Kaldunski, Anne E Kwitek, Stanley J F Laulederkind, Marek A Tutaj, 
  Mahima Vedi, Shur-Jen Wang, Peter D’Eustachio, Lucila Aimo, Kristian Axelsen, Alan Bridge, 
  Nevila Hyka-Nouspikel, Anne Morgat, Suzi A Aleksander, J Michael Cherry, Stacia R Engel, 
  Kalpana Karra, Stuart R Miyasato, Robert S Nash, Marek S Skrzypek, Shuai Weng, Edith D Wong, 
  Erika Bakker, Tanya Z Berardini, Leonore Reiser, Andrea Auchincloss, Kristian Axelsen, 
  Ghislaine Argoud-Puy, Marie-Claude Blatter, Emmanuel Boutet, Lionel Breuza, Alan Bridge, 
  Cristina Casals-Casas, Elisabeth Coudert, Anne Estreicher, Maria Livia Famiglietti, 
  Marc Feuermann, Arnaud Gos, Nadine Gruaz-Gumowski, Chantal Hulo, Nevila Hyka-Nouspikel, 
  Florence Jungo, Philippe Le Mercier, Damien Lieberherr, Patrick Masson, Anne Morgat, Ivo 
  Pedruzzi, Lucille Pourcel, Sylvain Poux, Catherine Rivoire, Shyamala Sundaram, Alex Bateman, 
  Emily Bowler-Barnett, Hema Bye-A-Jee, Paul Denny, Alexandr Ignatchenko, Rizwan Ishtiaq, Antonia 
  Lock, Yvonne Lussi, Michele Magrane, Maria J Martin, Sandra Orchard, Pedro Raposo, Elena 
  Speretta, Nidhi Tyagi, Kate Warner, Rossana Zaru, Alexander D Diehl, Raymond Lee, Juancarlos 
  Chan, Stavros Diamantakis, Daniela Raciti, Magdalena Zarowiecki, Malcolm Fisher, Christina 
  James-Zorn, Virgilio Ponferrada, Aaron Zorn, Sridhar Ramachandran, Leyla Ruzicka, Monte 
  Westerfield, The Gene Ontology knowledgebase in 2023, Genetics, Volume 224, Issue 1, May 2023, 
  iyad031, https://doi.org/10.1093/genetics/iyad031 
  - Ashburner, M., Ball, C., Blake, J. et al. Gene Ontology: tool for the unification of biology. 
  Nat Genet 25, 25–29 (2000). https://doi.org/10.1038/75556