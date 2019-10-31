import pandas as pd 
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


def calculateShannon(host):
    """
    Function for calculating Shannon Diveristy Index
    host: List 
    return int
    """
    sum = 0.0
    for species_abundance in host:
        if species_abundance != 0.0:
            sum += species_abundance * np.log(species_abundance)
    shannon_index = -sum 
    return shannon_index



data = pd.read_csv('MotherInfant.run_to_sample.txt', delim_whitespace = True)

run_acession_list = data['run_accession'].tolist()
cohort_list = data['cohort'].tolist()
cohorts = dict()
for index, run_accession in enumerate(run_acession_list):
    cohorts[run_accession] = cohort_list[index]

host_list = list(cohorts.keys())

species_data = pd.read_csv('relative_abundance.txt', delim_whitespace= True)
species_data = species_data.set_index("species_id", drop = True)

for host in list(species_data):
    if host not in host_list:
        del species_data[host]

host_shannon = dict()

for host in list(species_data):
    shannon_index = calculateShannon(species_data[host])
    host_shannon[host] = []
    host_shannon[host].append(shannon_index)
for host in cohorts:
    host_shannon[host].append(cohorts[host])

host_shannon_data = pd.DataFrame.from_dict(host_shannon, columns = ['Shannon Index', 'cohort'], orient='index')

w = 6
h = 6
d = 100
plt.figure(figsize=(w, h), dpi=d)
sns.boxplot(x="cohort", y="Shannon Index", data=host_shannon_data, order=["M", "B", "4M", "12M"])
sns.stripplot(x="cohort", y="Shannon Index", data=host_shannon_data, jitter=True, order=["M", "B", "4M", "12M"])
plt.show()
plt.savefig('mother_infant_species.pdf')  