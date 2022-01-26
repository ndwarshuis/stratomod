import pickle 

from dash import html
from interpret import set_visualize_provider
from interpret.provider import InlineProvider

import pandas as pd
from sklearn.model_selection import train_test_split

from interpret.glassbox import ExplainableBoostingClassifier
from interpret import show

df = pd.read_csv("HG002_DV_Illumina_df_ebm_indels_annotated_with_segdups_and_homopolymers_DP_VAF_indel_length_CG_AT_lens_label.bed", sep="\t")
df.columns = ["DP", "VAF", "indel_length", "CGhomopolgt3_len", "AThomopolgt3_len", "alignL_max", "alignL_count", "fracMatchIndel_max", "label"]

# adapted from #3 at https://towardsdatascience.com/converting-data-to-a-numeric-type-in-pandas-db9415caab0b
df['DP'] = pd.to_numeric(df['DP'], errors='coerce').fillna(0).astype('int')
df['VAF'] = pd.to_numeric(df['VAF'], errors='coerce').fillna(0.0).astype('float')
df['indel_length'] = pd.to_numeric(df['indel_length'], errors='coerce').fillna(0).astype('int')
df['CGhomopolgt3_len'] = pd.to_numeric(df['CGhomopolgt3_len'], errors='coerce').fillna(0).astype('int')
df['AThomopolgt3_len'] = pd.to_numeric(df['AThomopolgt3_len'], errors='coerce').fillna(0).astype('int')
df['alignL_max'] = pd.to_numeric(df['alignL_max'], errors='coerce').fillna(0).astype('int')
df['alignL_count'] = pd.to_numeric(df['alignL_count'], errors='coerce').fillna(0).astype('int')
df['fracMatchIndel_max'] = pd.to_numeric(df['fracMatchIndel_max'], errors='coerce').fillna(0.0).astype('float')

train_cols = df.columns[0:8]
label = df.columns[-1]
X = df[train_cols]
y = df[label]

seed = 1

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.20, random_state=seed)

ebm = ExplainableBoostingClassifier(random_state=seed)
ebm.fit(X_train, y_train)

with open('ebm_indels_data.pickle', 'wb') as f:
    pickle.dump(ebm, f, pickle.HIGHEST_PROTOCOL)

with open('indels_X_train_data.pickle', 'wb') as f:
    pickle.dump(X_train, f, pickle.HIGHEST_PROTOCOL)

with open('indels_X_test_data.pickle', 'wb') as f:
    pickle.dump(X_test, f, pickle.HIGHEST_PROTOCOL)

with open('indels_y_train_data.pickle', 'wb') as f:
    pickle.dump(y_train, f, pickle.HIGHEST_PROTOCOL)

with open('indels_y_test_data.pickle', 'wb') as f:
    pickle.dump(y_test, f, pickle.HIGHEST_PROTOCOL)