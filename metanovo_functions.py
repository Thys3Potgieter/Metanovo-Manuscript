#!/usr/bin/env python4

from functools import reduce
import pandas as pd
import scipy
import pandas as pd
import numpy as np
import os 
import shutil
import scipy.stats
import scikit_posthocs as ph
import numpy as np
import Bio
import scipy.stats
import scikit_posthocs as ph
import numpy as np
import matplotlib_venn
from matplotlib_venn import venn3, venn2
import matplotlib.pyplot as plt
from pylab import *
import matplotlib.pyplot as plt
from pylab import *
from IPython.display import Image
import skbio
import pandas as pd
import Bio
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from collections import defaultdict
from collections import Counter
import os
import matplotlib.pyplot as plt
from pylab import *
import pandas as pd
from IPython.display import Image
from IPython.display import HTML
from matplotlib_venn import venn3, venn2
from natsort import realsorted
import operator
from itertools import cycle, islice

import matplotlib
font = {'size'   : 16}

matplotlib.rc('font', **font)



def compare_msms(mm1, mm2, names, outpath):
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    print( "\nall peptide sequences")
    set1 = set(mm1['Sequence'].tolist())
    set2 = set(mm2['Sequence'].tolist())
    #figure(num=None, figsize=(8, 6), dpi=120, facecolor='w', edgecolor='k')
    fig = plt.figure()
    venn2([set1, set2], tuple(names))
    #plt.savefig(outpath + '/proteingroups_2pep.png')
    fig.savefig(outpath + '/all_peptide_sequences.png', bbox_inches='tight')
    plt.show()
    
    print("\nall peptide scan overlap")
    mset1 = set(mm1['ScanID'].tolist())
    mset2 = set(mm2['ScanID'].tolist())
    #figure(num=None, figsize=(8, 6), dpi=120, facecolor='w', edgecolor='k')
    fig = plt.figure()
    venn2([mset1, mset2], tuple(names))
    fig.savefig(outpath + '/all_peptide_scans.png', bbox_inches='tight')

    #fig.savefig('analysis/mli_venn1.png', bbox_inches='tight')
    plt.show()
    
    print('\nshared peptide scan overlap')
    exclusive_set1 = set1 & set2
    m1  = set(mm1[mm1['Sequence'].isin(exclusive_set1)]['ScanID'].tolist())
    m2  = set(mm2[mm2['Sequence'].isin(exclusive_set1)]['ScanID'].tolist())
    #figure(num=None, figsize=(8, 6), dpi=120, facecolor='w', edgecolor='k')
    fig = plt.figure()
    venn2([m1, m2], tuple(names))
    #fig.savefig('analysis/mli_venn1.png', bbox_inches='tight')
    fig.savefig(outpath + '/shared_peptide_scans.png', bbox_inches='tight')

    plt.show()
    
    print('\nexclusive peptide scan overlap')
    exclusive_set1 = set1 - set2
    m1  = set(mm1[mm1['Sequence'].isin(exclusive_set1)]['ScanID'].tolist())
    
    exclusive_set2 = set2 - set1
    m2  = set(mm2[mm2['Sequence'].isin(exclusive_set2)]['ScanID'].tolist())
    #figure(num=None, figsize=(8, 6), dpi=120, facecolor='w', edgecolor='k')
    fig = plt.figure()
    venn2([m1, m2], tuple(names))
    fig.savefig(outpath + '/exclusive_peptide_scans.png', bbox_inches='tight')

    #fig.savefig('analysis/mli_venn1.png', bbox_inches='tight')
    plt.show()

def process_txt(path, name, outfolder, col_order = None, prot_map = None, dropna=False):
    print('Analysis: {}'.format(name))
    analysis_folder = outfolder
    if not os.path.exists(analysis_folder):
        os.mkdir(analysis_folder)
    print(path)
    result = {}
    peptides = pd.read_csv('{}/peptides.txt'.format(path), sep='\t', low_memory=False)
    proteins = pd.read_csv('{}/proteinGroups.txt'.format(path), sep='\t', low_memory=False)
    
    try:
        summary = pd.read_csv('{}/summary.txt'.format(path), sep='\t', low_memory=False)
    except:
        pass
    target_peptides = peptides[peptides['Reverse']!= '+']
    reverse_peptides = peptides[peptides['Reverse'] == '+']
    target_proteins = proteins[proteins['Reverse']!= '+']
    
    if not prot_map is None:
        target_proteins = proteingroup_organisms(target_proteins, prot_map, dropna = dropna)
        target_peptides = peptide_organisms(target_peptides, target_proteins, dropna = dropna )
    
    print('{}: Total peptides: '.format(name),len(peptides))
    print('{}: Total target peptides: '.format(name),len(target_peptides))
    print('{}: Total proteins: '.format(name),len(proteins))
    print('{}: Total target proteins: '.format(name),len(target_proteins))
    try:
        print('{}: Total MS/MS Identified [%]: '.format(name), summary[-1:]['MS/MS Identified [%]'].values[0])
    except:
        pass
    peptide_sequences = set(target_peptides['Sequence'].tolist())
    
    if os.path.exists('{}/evidence.txt'.format(path)):
        evidences = pd.read_csv('{}/evidence.txt'.format(path), sep='\t', low_memory=False)
        counts = spectral_counts(evidences, target_peptides, col_order)
        result['SpectralCounts'] = counts
    
    result['TargetPeptides'] = target_peptides
    result['ReversePeptides'] = reverse_peptides
    result['TargetProteins'] = target_proteins
    result['Path'] = path

    reverse_peptide_sequences = reverse_peptides['Sequence'].tolist()
    w = open('{}/peptide_list.txt'.format(analysis_folder),'w')
    w.write('\n'.join(peptide_sequences))
    w.close()
    
    return result

def list_kw_dunn(names, data, value, group, path):
    colnames=names
    kw = scipy.stats.kruskal(*data)
    w = open(path + '/kw.txt', 'w')
    w.write(str(kw))
    print(kw)
    w.close()
    post_hoc = pd.DataFrame(ph.posthoc_dunn(data, p_adjust = 'fdr_bh'))
    post_hoc.index = names
    post_hoc.columns = names
    post_hoc.to_csv(path + '/dunn_bh.csv')
    return post_hoc

def get_organism(rec):
    org = rec.description.split('OS=')[1]
    if '=' in org:
        org = org.split('=')[0].split()[:-1]
        org = ' '.join(org)
    return org

def get_protein(rec):
    prot = rec.id.split()[0]
    return prot
    
def get_organisms(fasta):
    orgs = defaultdict(int)
    for rec in fasta:
        org = get_organism(rec)
        orgs[org] += 1
    return orgs

def prot2organism(fasta):
    prots = {}
    for rec in fasta:
        org = get_organism(rec)
        prot = get_protein(rec)
        prots[prot] = org
    return prots

def create_fasta(peptide_list):
    count = 1
    recs = []
    for pep in peptide_list:
        rec = SeqRecord(id = 'peptide_{}'.format(str(count)), seq = Seq(pep))
        count += 1
        recs.append(rec)
    return recs

def id2organisms(val,  map):
    orgs = []
    for id in val.split(';'):
        orgs.append(map[id])
    orgs = list(set(orgs))
    orgs.sort()
    return ';'.join(orgs)

def proteingroup_organisms(pg, protein_mapping, dropna=False):
    pg['LeadingProtein'] = pg['Protein IDs'].apply(lambda x : x.split(';')[0])
    pg['OS'] = pg['LeadingProtein'].map(protein_mapping)
    if dropna == True:
        pg = pg[pg['OS'].notnull()]
    return pg

def peptide_organisms(peptides, proteingroups, dropna=False, unique=True):
    #if unique == True
    #peptides = peptides[peptides['Protein group IDs'].apply(lambda x : x.find(';') == -1)]
    proteingroups['id'] = proteingroups['id'].apply(str)
    pepmap = proteingroups.set_index('id')['OS']
    pepmap_prot = proteingroups.set_index('id')['LeadingProtein']
    peptides['OS'] = peptides['Protein group IDs'].map(pepmap)
    peptides['LeadingProtein'] = peptides['Protein group IDs'].map(pepmap_prot)
    if dropna == True:
        peptides = peptides[peptides['OS'].notnull()]
    return peptides

def evidence_organisms(evidences, peptides, dropna=False):
    pepmap = peptides.set_index('Sequence')['OS']
    evidences['OS'] = evidences['Sequence'].map(pepmap)

    pepmap_prot = peptides.set_index('Sequence')['LeadingProtein']
    evidences['LeadingProtein'] = evidences['Sequence'].map(pepmap_prot)

    if dropna == True:
        evidences = evidences[evidences['OS'].notnull()]
    return evidences


def spectral_counts(evidences, peptides, order):
    counts = evidences.groupby(['Sequence', 'Experiment']).agg({"MS/MS Count":'sum'})
    counts = counts.reset_index('Experiment')
    counts['Experiment' ] = counts['Experiment'].apply(lambda x : 'Counts {}'.format(x))
    count_cols = order
    counts = counts.pivot(columns='Experiment', values = 'MS/MS Count').fillna(0)
    counts = counts.reset_index('Sequence')
    if not count_cols is None:
        newcounts = counts[count_cols]
        newcounts['Sequence'] = counts['Sequence']
    else:
        newcounts = counts
    if 'OS' in peptides.columns:
        pepmap = peptides.set_index('Sequence')['OS']
        newcounts['OS'] =  newcounts['Sequence'].map(pepmap)
    if 'LeadingProtein' in peptides.columns:
        pepmap = peptides.set_index('Sequence')['LeadingProtein']
        newcounts['LeadingProtein'] =  newcounts['Sequence'].map(pepmap)

    newcounts = newcounts[newcounts['Sequence'].isin(peptides['Sequence'].tolist())]
    return newcounts


def get_top_match_blast(table_path):
    df = pd.read_csv(table_path)
    df = df.sort_values('hsp.score', ascending=False)
    df = df.drop_duplicates('blast_record.query', keep='first')
    print(df.columns.tolist())
    return df
    

def plot_dataframe_bar(df, title, outpath, legend=True):
    f = plt.figure()
    plt.title(title, color='black')
    df.plot(kind='bar', ax=f.gca(), legend=legend)
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    if legend == False:
        plt.gca().legend_.remove()

    plt.show()
    plt.savefig(outpath)

def organism_pep(peptides, organism_map):
    df = peptides[peptides['Unique (Groups)'] == 'yes' ]
    total = sum(list(organism_map.values()))
    probs = {}
    for org, table in df.groupby('OS'):
        org_counts = organism_map[org]
        org_prior = org_counts / total
        peps = table['PEP'].tolist()
        prob = reduce( operator.mul , peps, 1) * org_prior
        probs[org] = prob
    probs = pd.DataFrame.from_dict(probs, orient='index')
    probs.rename(columns={0:'PEP','0':'PEP'}, inplace=True)
    probs = probs.sort_values('PEP', ascending=True)
    return probs

def analyze_organisms(pg , peptides , counts, organism_counts, outpath):
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    df = pg
    print("Protein group analysis:\n")
    print("Number of protein groups: {}".format(len(pg)))
    df['LeadingProteins'] = 1
    olp = df.groupby('OS').agg({'LeadingProteins' : 'sum'}).rename(columns={"LeadingProteins":"Leading proteins with 1 or more peptides"})
    print(olp)

    df2 = df[df['Unique peptides'] >= 2]
    upc2 = df2.groupby('OS').agg({'LeadingProteins' : 'sum'}).rename(columns={"LeadingProteins":"Leading proteins with 2 or more peptides"})
    print(upc2)
    
    #upc = df.groupby('OS').agg({'Score' : 'median'}).rename(columns={"Score":"Median protein group score"})
    #print(upc)
    #upc.plot.bar()
    
    up = df.groupby('OS').agg({'Unique peptides' : 'sum'})
    print(up)
    
    mup = df.groupby('OS').agg({'Unique peptides' : 'mean'}).rename(columns={"Unique peptides":"Mean protein group unique peptide count"})
    print(mup)
   
    del df['LeadingProteins']
    
    # Peptides
    print("\nPeptide analysis:\n")
    print("Number of peptides: {}".format(len(peptides)))

    peptides['SequenceCount'] = 1
    olp_pep = peptides.groupby('OS').agg({'SequenceCount' : 'sum'}).rename(columns={"SequenceCount":"Peptide count"})
    print(olp_pep)
   
    upc_pep = peptides.groupby('OS').agg({'PEP' : 'median'}).rename(columns={"PEP":"Median peptide PEP score"})
    #print(upc_pep)
    print()
    
    # Organism PEP
    o_pep = organism_pep(peptides, organism_counts)
    
    f= outpath + '/organism_pep.csv'
    o_pep.to_csv(f)
    print(f)
    print(o_pep)
    print()

    # Spectral Counts
    print("\nSpectral count analysis:\n")
    agg_dict = {}
    count_cols = []
    for col in counts.columns:
        if col.startswith('Counts '):
            count_cols.append(col)
            agg_dict[col] = 'sum'
    pep_msms = counts.groupby('OS').agg(agg_dict)[count_cols]
    org_msms = pep_msms.transpose()
   
    upc_pep = peptides.groupby('OS').agg({'PEP' : 'median'}).rename(columns={"PEP":"Median peptide PEP score"})
    #print(upc_pep)
    print()
    
    #f = outpath + '/organism_pep.png'
    #plot_dataframe_bar(o_pep , 'Organism PEP scores', f)

    f = outpath + '/sample_counts.png'
    plot_dataframe_bar(org_msms, 'Organism spectral counts', f)

    f = outpath + '/proteingroups.png'
    plot_dataframe_bar(olp, 'Leading proteins with one or more peptides', f, False)
    
    f = outpath + '/proteingroups_2pep.png'
    plot_dataframe_bar(upc2, 'Leading proteins with two or more peptides', f, False)
    
    f = outpath + '/unique_peptides.png'
    plot_dataframe_bar(up, 'Protein group unique peptide counts', f, False)
        
    f = outpath + '/mean_unique_peptides.png'
    plot_dataframe_bar(mup, 'Mean protein group unique peptide counts', f, False)
    
    f = outpath + '/peptide_count.png'
    plot_dataframe_bar(olp, 'Peptide count', f, False)
        

    pep_figure(peptides, outpath + '/organism_pep.png', 'OS')
    
    #up = df.groupby('OS').agg({'Unique peptides' : 'sum'})
    #print(up)
    del peptides['SequenceCount']

    


def pep_figure(peptides, outpath, col):
    #figure(num=None, figsize=(16, 12), dpi=120, facecolor='w', edgecolor='k')
    names = []
    pep_scores = []
    for name, group in  peptides.groupby(col):
        names.append(name)
        pep_scores.append(group['PEP'].tolist())
        
    colours = ['b','g','r','c','m','y','k', 'b','g','r','c','m','y','k', 'b','g','r','c','m','y','k']
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Create the boxplot
    bp = ax.boxplot(pep_scores, patch_artist=True, showfliers=False)
    ## change outline color, fill color and linewidth of the boxes
    count = 0
    col_ind=0
    for box in bp['boxes']:
        count += 1
        # change outline color
        box.set( color='#7570b3', linewidth=1)
        # change fill color
        box.set( facecolor = colours[col_ind] )
        ## change color and linewidth of the whiskers
        if count % 3 == 0:
            col_ind +=1
    for whisker in bp['whiskers']:
        whisker.set(color='#7570b3', linewidth=1)
    ## change color and linewidth of the caps
    for cap in bp['caps']:
        cap.set(color='#7570b3', linewidth=1)
    ## change color and linewidth of the medians
    for median in bp['medians']:
        median.set(color='#b2df8a', linewidth=2)
        #median.set(linewidth=2)
    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='.', color='#e7298a', alpha=0.5)
    ## Custom x-axis labels
    ax.set_xticklabels(names, rotation=90)
    #ax.set_yticklabels('Posterior Error Probability (PEP)') 
    ax.set_title('Peptide PEP Score distributions')
    ## Remove top axes and right axes ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    fig.savefig( outpath , bbox_inches='tight')
    plt.show()
    fig.clf()
    return
  
def peps2blast(peptides, blast):
    _ = peptides.copy()
    _ = _[['Sequence','OS']].sort_values('OS').drop_duplicates('Sequence')
    _ = _.rename(columns = {'OS':"ProteinGroup_OS"})
    _ = _.reset_index()
    del _['index']
    _ = pd.merge(_, blast, how='outer', left_on='Sequence', right_on='_query.sequence')
    return _

def get_scannum(msms):
    #msms = msms[msms['Sequence'].isin(peptide_list)]
    msms['ScanID'] = msms['Raw file'].apply(str) + '_' + msms['Scan number'].apply(str) + '_' + msms['Scan index'].apply(str) 
    return msms 

def plot_taxa_runs( data_results, col , run, mapping, how='left'):
    #keys = list(data_results.keys())
    plt.clf()
    # Peptide counts per sample
    count_df = pd.DataFrame()
    target_peptides = data_results[run]['SpectralCounts']
    dfs = []
    for c in target_peptides.columns:
            if c.startswith('Counts '):
                seqs = target_peptides[[ 'Sequence' , c ]]
                seqs = seqs[seqs[c] > 0]
                mapped = pd.merge(seqs, mapping, how=how, left_on='Sequence', right_on='peptide')
                agg_cols = { c: sum }
                agg = mapped.groupby(mapped[col]).agg(agg_cols)
                #print(agg)
                agg[c] = agg[c].apply(int)
                agg.rename(columns={c:c.split('Counts ')[1]},inplace=True)
                dfs.append(agg)
    count_df = pd.concat(dfs, axis=1)            
    count_df = count_df.replace(np.nan, 0)
    count_df = count_df.astype(int64)
    count_df['Total'] = count_df.sum(axis=1)
    count_df = count_df.sort_values('Total',ascending=False)
    del count_df['Total']
    count_df_all = cound_df.copy()
    count_df = count_df[:15]
    t_df = count_df.transpose()    
    ax1 = t_df.plot(kind='line', rot=1 )#, stacked=True)
    ax1.set_title("Spectral Count by UniPept pept2lca {}".format(col))
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
              ncol=3, fancybox=True, shadow=True)
    fig = ax1.get_figure()
    plt.xticks(rotation=90)
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.show()
    #fig.savefig('analysis/{}_bar.png'.format(col), bbox_inches='tight')
    print(len(count_df_all))
    return count_df_all

def plot_taxa( results, col , mapping, order, filt=0.5, taxon_filter={}, truth_list=None , level='msms', truth_exclude=False, filter_method='include', fname=None): 
    keys = list(results.keys())
    plt.clf()
    # Peptide counts per sample
    count_df = pd.DataFrame()
    mapping_copy = mapping.copy()
    
    # exclude peptides not quantified at a given level
    mapping_copy = mapping_copy[mapping_copy[col].notnull()]
    if not truth_list is None:
        rename = {}
        possible = set(mapping_copy[col].dropna().tolist())
        for _ in possible:
            if not _ in truth_list:
                if truth_exclude == False:
                    rename[_] = 'missasignment'
                else:
                    rename[_] = _
            else:
                if truth_exclude==False:
                    rename[_] = _
                else:
                    rename[_] = 'correct assignment'
        

    if level =='msms':
        agglevel = 'MS/MS Count'
    elif level == 'peptide':
        agglevel = 'Peptide Count'

    for name in order:
        target_peptides = results[name]['TargetPeptides']
        mapped = pd.merge(target_peptides, mapping_copy, how='inner', left_on='Sequence', right_on='peptide') 
        mapped['Peptide Count'] = 1 # so we can aggregate peptide counts per taxon level
        #print(mapped.head())
        #mapped = mapped.replace(np.nan, 'Uncharacterized')
        agg_cols = {'MS/MS Count': sum, 'Peptide Count': sum}
        #print('Taxon filter: ', taxon_filter)
        for _ in taxon_filter:
            agg_cols[_] = 'first'
        agg = mapped.groupby(mapped[col]).agg(agg_cols)
        #print("agg cols: ", agg_cols)
        #print('Agg: ', agg) 
        agg['Peptide Count %'] = agg['Peptide Count']  / agg['Peptide Count'].sum() * 100.0
        
        #print(set(agg['superkingdom_name'].tolist()))
        agg = agg[agg['Peptide Count %'] >= filt]
        #print(set(agg['superkingdom_name'].tolist()))
        
        for _ in taxon_filter:
            #print(_)
            if filter_method=='include':
                agg = agg[agg[_] == taxon_filter[_]]
            elif filter_method=='exclude':
                agg = agg[agg[_] != taxon_filter[_]]
        count_df[name ] = agg[agglevel]
    #print(count_df)
    agg_cols = {}
    for _ in count_df.columns:
        agg_cols[_] =  sum
    count_df = count_df.reset_index()
    if truth_list is not None:
        count_df[col] = count_df[col].map(rename)
    count_df = count_df[count_df[col] != 'correct assignment']
    #print(count_df) 
    count_df = count_df.groupby(count_df[col]).agg(agg_cols)
    #print(agg)
    #count_df = count_df.set_index(col)
    
    count_df = count_df.replace(np.nan, 0)
    count_df = count_df.astype(int)
    count_df['Max'] = count_df.max(axis=1)
    count_df = count_df.sort_values('Max',ascending=False)
    del count_df['Max']
    count_df_all = count_df.copy()
    limit=15
    other = count_df[limit:].sum()
    
  
    count_df = count_df[:limit]
    if sum(other.values) > 0:
        count_df.loc['other'] = other
    
    for column in count_df_all.columns:
        vals = count_df_all[count_df_all[column] > 0 ]
        #print(column, len(vals))
   
    my_colors = list(islice(cycle(['skyblue','plum','navajowhite','hotpink', 'lavender','brown', 'limegreen', 'green', 'darkorchid','papayawhip','mediumblue','slategrey','tan','indigo', 'thistle','teal']), None, len(count_df)))
    ind = list(count_df.index)
    pos=-1
    if 'missasignment' in ind:
        ind = [i for i in ind if not i == 'missasignment'] + ['missasignment']

        count_df = count_df.loc[ind]

        pos = ind.index('missasignment')
        my_colors[pos] = 'r'
    t_df = count_df.transpose()    
    #matplotlib.rcParams.update({'font.size': 11})
    ax1 = t_df.plot(kind='bar', rot=1, stacked=True, color=my_colors,figsize=(10,6)) 
    _ = {}
    _count = {}
    _total = {}
    for p in ax1.patches:
        x = p.get_x()
        if not x in _total:
            _total[x] = 0
        _total[x] += p.get_height()

    for p in ax1.patches:
        #print(dir(p))
        #print(p.get_name())
        x = p.get_x()
        if not x in _:
            _[x] = 0
            _count[x] = 0
        #if _count[x] == pos:
        height = p.get_height()
        
        if _total[x] == 0:
            continue
        pct = height / _total[x]  * 100.0
        #if pct < 1:
        #    continue
        pct = np.round(pct, 2)
        #if pct > 1:
        #    ax1.annotate(str(height) + ' (' + str(pct) + ' %)', (p.get_x() * 1.005, _[x] * 1.005))
        _[x] += p.get_height()
        _count[x] += 1
    if truth_exclude == True:
        heading = 'False taxa {}'.format(agglevel)
    else:
        heading = agglevel
    heading = heading[:1].upper() + heading[1:].lower()
    ax1.set_title("{} by UniPept pept2lca {}".format(heading, col))
    ax1.legend(bbox_to_anchor=(1.05, 0), loc='lower left', borderaxespad=0.)
    plt.xticks(rotation=45)
    plt.xlabel('Database')
    plt.ylabel('{} count'.format(level[:1].upper() + level[1:]))
    fig = ax1.get_figure()
    if not fname is None:
        fig.savefig(fname, bbox_inches='tight', dpi=600)
    plt.show()
    #print(count_df)
    #count_df_all.loc['Total'] = count_df_all.sum()
    return count_df_all


from PIL import Image, ImageDraw, ImageFont
def concat_4way_image(tl, tr, bl, br, outfile):
    fnt  = ImageFont.truetype("/cbio/users/ptgmat003/2018/Metaproteomics/Thys/metanovo/Arial.ttf", 150, encoding="unic")

    im_tl = Image.open(tl)
    d = ImageDraw.Draw(im_tl)
    d.multiline_text((0,0), "A", font=fnt, fill=(0, 0, 0))

    im_tr = Image.open(tr)
    d = ImageDraw.Draw(im_tr)
    d.multiline_text((0,0), "B", font=fnt, fill=(0, 0, 0))

    im_bl = Image.open(bl)
    d = ImageDraw.Draw(im_bl)
    d.multiline_text((0,0), "C", font=fnt, fill=(0, 0, 0))
   
    im_br = Image.open(br)
    d = ImageDraw.Draw(im_br)
    d.multiline_text((0,0), "D", font=fnt, fill=(0, 0, 0))
    
    mid_x = max(im_bl.width, im_tl.width)
    mid_y = max(im_tl.height, im_tr.height)

    top = Image.new('RGB', (mid_x + im_tr.width, mid_y), color=(255,255,255))
    top.paste(im_tl, (0, 0))
    top.paste(im_tr, (mid_x, 0))
    
    bottom = Image.new('RGB', (mid_x + im_br.width, max(im_br.height, im_bl.height)), color=(255,255,255))
    bottom.paste(im_bl, (0, 0))
    bottom.paste(im_br, (mid_x, 0))

    dst = Image.new('RGB', (max(top.width, bottom.width), top.height + bottom.height), color=(255, 255, 255))
    dst.paste(top, (0, 0))
    dst.paste(bottom, (0, top.height))
    dst.save(outfile)
    return dst


