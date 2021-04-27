import pandas as pd


HEADER = """
[Header],,,,,,,,,
IEMFileVersion,1,,,,,,,,
Investigator Name,Alexander Misharin,,,,,,,,
Experiment Name,Aging_PF_ISRIB,,,,,,,,
Date,,,,,,,,,
Workflow,GenerateFASTQ,,,,,,,,
Application,NextSeq FASTQ Only,,,,,,,,
Assay,TruSeq HT,,,,,,,,
Description,,,,,,,,,
Chemistry,Amplicon,,,,,,,,
,,,,,,,,,
[Reads],,,,,,,,,
,,,,,,,,,
,,,,,,,,,
,,,,,,,,,
[Settings],,,,,,,,,
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,,,,,,,,
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT,,,,,,,,
[Data],,,,,,,,,
SampleID,Sample_Name,Sample_Plate,Organism,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
"""


def reverse_tr(x):
    return x[::-1].translate(str.maketrans({'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A'}))


def prepare_sample_sheet(samples_file, output, project):
    columns = [
        "SampleID",
        "Name",
        "Plate",
        "Species",
        "Index1Name",
        "Index1Sequence",
        "Index2Name",
        "Index2Sequence",
        "Project",
        "NucleicAcid",
    ]

    samples = pd.read_table(samples_file)
    samples.Species = "mm10"
    samples.Project = project
    samples.Index2Sequence = samples.Index2Sequence.apply(reverse_tr)
    samples.NucleicAcid = samples.NucleicAcid + "-seq"
    samples["Plate"] = ""

    with open(output, "w") as f:
        f.write(HEADER.lstrip())
        f.write(samples.loc[:, columns].to_csv(None, index=False, header=False))
