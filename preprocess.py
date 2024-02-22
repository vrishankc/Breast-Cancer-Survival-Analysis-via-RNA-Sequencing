import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--clinical', type=str,
                    help='path to the clinical data csv')
parser.add_argument('--rna', type=str,
                    help='path to rna csv data')
parser.add_argument('--pam', type=str, help = 'path to PAM50 protein assay data')




## PREPROCESSING FOR RNA

def extract_prefix_with_tag(phrase_list):
    prefixes_with_tag = []
    for phrase in phrase_list:
        prefix = phrase.split('.')[0]
        prefix_with_tag = "TCGA-" + prefix.replace("–", "-")
        prefixes_with_tag.append(prefix_with_tag)
    return prefixes_with_tag

if __name__ == '__main__':
  clinical = pd.read_csv(args.clinical)
  rna = pd.read_csv(args.rna)
  pam = pd.read_csv(args.pam)

  new_cols = extract_prefix_with_tag(rna.columns.tolist())
  new_cols[0] = "RefSeq_accession_number"
  new_cols[1] = "gene_symbol"
  new_cols[2] = "gene_name"
  rna.columns = new_cols
  rna = rna.transpose()
  gene_refseq = rna.iloc[0:2]
  rna = rna.drop(rna.index[[1, 2]], axis=0)
  rna = rna.transpose().dropna()
  rna = rna.transpose()
  rna.head(10)
  rna_ids = rna.transpose().columns[1:].tolist()

  ## PREPROCESSING FOR CLINICAL
  clinical = clinical[["Complete TCGA ID", "Gender", "Age at Initial Pathologic Diagnosis", "OS Time", "OS event"]]
  clinical.dropna()
  
  for ids in rna_ids:
      if ids in clinical["Complete TCGA ID"].tolist():
          time = clinical.loc[clinical["Complete TCGA ID"] == ids, "OS Time"].values
          event = clinical.loc[clinical["Complete TCGA ID"] == ids, "OS event"].values
          rna.loc[rna.index == ids, "time"] = time[0] if len(time) > 0 else None
          rna.loc[rna.index == ids, "event"] = event[0] if len(event) > 0 else None
      else:
          print(ids)
  rna = rna.dropna()

  rna.to_csv("rna.csv")
