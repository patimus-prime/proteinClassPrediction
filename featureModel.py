import pandas as pd
import sys as sys
import os as os
import ast as ast
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem import PandasTools
from rdkit import RDConfig
from rdkit.Chem import PandasTools as pdt
from rdkit.Chem import Descriptors as D

big_df = pd.read_csv("smiles_l1.csv")  # warning: BIG file

# 1 million (Dr.Evil moment) samples! let's get 1%, check still get 15 types
df = big_df.sample(frac=0.01, replace=False, random_state=1)
df.head()

# get those lists in the 2nd column resolved; unknown to me still if they're just dup'd
df.L1_class_name = df.L1_class_name.apply(ast.literal_eval)

# df.L1_class_name.unique()
# df['L1_class_name'].unique().sort()

# now 1 element lists - still need to convert to a string so we can factorize
df['string_L1'] = df['L1_class_name'].apply(lambda x: x[0])
df.head()
# so this lambda fn; we set x = the first and only element in L1_class_name; that
# x becomes string_L1 in our new column. lambda will make a continued appearance
# during RDKit feature generation

# pd DEMANDS you declare your factors and get the KEY of what they are
df['Target_Class'], uniques = pd.factorize(df.string_L1)
df.head()
#print("V important, the uniques!!!: ")
# uniques

pdt.AddMoleculeColumnToFrame(df, 'canonical_smiles', includeFingerprints=True)

# get some descriptions of the data, this one is SA:
# D.TPSA(df.ROMol.iat[0])
# D.MolLogP(df.ROMol.iat[0])
# D.MolWt(df.ROMol.iat[0])


df['TPSA'] = df['ROMol'].apply(lambda x: D.TPSA(x))
df['MolLogP'] = df['ROMol'].apply(lambda x: D.MolLogP(x))
df['MolWt'] = df['ROMol'].apply(lambda x: D.MolWt(x))

# getting ahead of ourselves, but declare the salt remover
# Cl, Br, Zn, Mg, Ca, Sr, Na, K
remover = SaltRemover(defnData="[Cl,Br,Zn,Mg,Ca,Sr,Na,K]")
