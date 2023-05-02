
import argparse as ap
import csv
import json
import numpy as np
import os
import pybedtools
import scipy

def loadSampleXciCpgs(prefix, xci_bed_fn):
    '''
    This will load all the XCI loci for a sample into memory for downstream analysis.
    :param prefix: - the prefix of the input files from tool `aligned_bam_to_cpg_scores`
    :param xci_bed_fn: - the BED file containing XCI loci
    :return: - a nested dictionary where the first key is the position on chrX and the second key is 'combined', 'hap1', or 'hap2'
        followed by data for that locus / source from the data load
    '''
    # first load the XCI bed file
    xci_bed = pybedtools.BedTool(xci_bed_fn)

    # now load the sample files
    file_loads = {
        # we do not need 'combined' key anymore
        # 'combined' : f'{prefix}.combined.bed.gz',
        'hap1' : f'{prefix}.hap1.bed.gz',
        'hap2' : f'{prefix}.hap2.bed.gz'
    }
    key_order = ['hap1', 'hap2']

    # return dictionary: ret[variant_position][key_type]
    # variant_position - a coordinate on chrX
    # key_type - one of the values in "key_order"
    ret = {}
    
    # now load each one
    for key_type in key_order:
        # file checks
        filename = file_loads[key_type]
        print(f'Loading {filename}...')
        if not os.path.exists(filename):
            if key_type in ['hap1', 'hap2']:
                raise Exception(f'File does not exist, was the BAM haplotagged? : {filename}')
            else:
                raise Exception(f'File does not exist, is prefix correct? : {filename}')
        
        # load the sample bed and intersect with the XCI coordinates
        sample_bed = pybedtools.BedTool(filename)
        print('\tIntersecting with XCI coordinates...')
        intersection = sample_bed.intersect(xci_bed)

        # iterate over the intersection with some sanity checks
        print('\tParsing intersection result...')
        for row in intersection:
            # this checks the field length to determine the pileup type
            if len(row.fields) == 10:
                # "count" pileup type, the official supported type
                chrom, start, stop, methyl_prob, haplotype, read_count, meth_count, unmeth_count, avg_meth_score, avg_unmeth_score = row.fields

                # dummy values
                discrete_methyl_prob = methyl_prob
            else:
                # "model" pileup type, technically works but we do not support
                chrom, start, stop, methyl_prob, haplotype, read_count, meth_count, unmeth_count, discrete_methyl_prob = row.fields

                # dummy values
                avg_meth_score = 1.0
                avg_unmeth_score = 0.0
            
            read_count = int(read_count)
            
            # sanity checks
            meth_count = int(meth_count)
            avg_meth_score = float(avg_meth_score)
            unmeth_count = int(unmeth_count)
            avg_unmeth_score = float(avg_unmeth_score)
            try:
                assert(chrom == 'chrX')
                assert(meth_count == 0 or avg_meth_score > 0.5)
                assert(unmeth_count == 0 or avg_unmeth_score < 0.5)
            except Exception as e:
                print(row.fields)
                raise e

            # things seem good, add to the list
            start = int(start)
            result_prob = float(methyl_prob) / 100.0
            
            if start not in ret:
                ret[start] = {}
            assert(key_type not in ret[start])

            # add to our dictionary by start position and key type
            ret[start][key_type] = {
                'read_count' : read_count,
                'methyl_count' : meth_count,
                'unmethyl_count' : unmeth_count,
                'methyl_prob': result_prob 
            }

    return ret

def calculateStatistics(sample_cpgs):
    '''
    This will parse the loaded CpG values and calculate various statistics on the results
    :param sample_cpgs: - the loaded sample XCI cpg values
    :return: - tuple (xci_stats, data_arrays)
        xci_stats - a dictionary of statistics related to XCI skew
        data_arrays - a dictionary of numpy arrays forming ordered lists of the loaded datasets, primarily for downstream visualizations
    '''
    # gather h1 and h2 into some arrays
    h1_vals = []
    h2_vals = []
    h1_readcount = []
    h2_readcount = []
    h1_unmethyl_reads = []
    h1_methyl_reads = []
    h2_unmethyl_reads = []
    h2_methyl_reads = []

    for k in sample_cpgs:
        # print(sample_cpgs[k])
        h1 = sample_cpgs[k].get('hap1', None)
        h2 = sample_cpgs[k].get('hap2', None)

        # make sure we have an entry for both haplotypes
        if h1 != None and h2 != None:
            h1_mp = h1['methyl_prob']
            h2_mp = h2['methyl_prob']

            h1_rcount = h1['read_count']
            h2_rcount = h2['read_count']

            h1_ucount = h1['unmethyl_count']
            h1_mcount = h1['methyl_count']
            h2_ucount = h2['unmethyl_count']
            h2_mcount = h2['methyl_count']

            # comparing methylation probability of the two haplotypes
            h1_vals.append(h1_mp)
            h2_vals.append(h2_mp)
            h1_readcount.append(h1_rcount)
            h2_readcount.append(h2_rcount)

            # pileup counts of both
            h1_unmethyl_reads.append(h1_ucount)
            h1_methyl_reads.append(h1_mcount)
            h2_unmethyl_reads.append(h2_ucount)
            h2_methyl_reads.append(h2_mcount)
    
    # convert to numpy arrays for manipulations
    h1_vals = np.array(h1_vals)
    h2_vals = np.array(h2_vals)

    # figure out the max and min of each h1/h2 pair
    hap_stack = np.vstack([h1_vals, h2_vals])
    argmin_hap = np.argmin(hap_stack, axis=0)
    argmax_hap = 1 - argmin_hap
    
    # build max/min haps
    min_hap = hap_stack[argmin_hap, np.arange(hap_stack.shape[1])]
    max_hap = hap_stack[argmax_hap, np.arange(hap_stack.shape[1])]
    
    # we now have ordered points such that we can do cluster stats
    cluster_points = np.vstack([min_hap, max_hap]).transpose()
    
    # weights are now assigned based on the actual observed min or max methylation prob
    rc_stack = np.vstack([h1_readcount, h2_readcount])
    min_weights = rc_stack[argmin_hap, np.arange(hap_stack.shape[1])]
    max_weights = rc_stack[argmax_hap, np.arange(hap_stack.shape[1])]
    
    # min_weights then max_weights to match construction of cluster_points
    rc_weights = np.vstack([min_weights, max_weights]).transpose()
    
    # now all the summary stats
    cluster_weighted_means = np.average(cluster_points, weights=rc_weights, axis=0)
    cluster_means = np.mean(cluster_points, axis=0)
    cluster_medians = np.median(cluster_points, axis=0)
    
    # statistical calculations
    combined_weight = np.sum(rc_weights, axis=1)
    delta_scores = max_hap - min_hap
    delta_mean = np.mean(delta_scores)
    weighted_mean = np.average(delta_scores, weights=combined_weight)
    delta_median = np.median(delta_scores)
    
    # lin regression mode:
    # logical form            ==>  standard linear alg. form for solving for p and e
    # frac + e = obs_low      ==>  frac  + e = obs_low
    # 1 - frac + e = obs_high ==>  -frac + e = obs_high - 1
    # this stores the constants in front of the unknown variable on the left-hand side
    systemLHS = [[1, 1],
                 [-1, 1]]
    
    # this stores the calculate constant on the right-hand side based on the observed average allelic ratios
    # systemRHS = [np.mean(min_hap), np.mean(max_hap) - 1]
    systemRHS = [np.median(min_hap), np.median(max_hap) - 1]
    
    # calculate and return frac, e
    linalg_result = np.linalg.solve(systemLHS, systemRHS)

    # we can also do this for all of our haps
    systemLHS = np.array([[1, 1]] * len(min_hap) + [[-1, 1]] * len(max_hap))
    systemRHS = np.hstack([min_hap, (max_hap - 1)])
    lstsq_result = np.linalg.lstsq(systemLHS, systemRHS, rcond=None)[0]

    # weighted form, via https://stackoverflow.com/questions/27128688/how-to-use-least-squares-with-weight-matrix
    flat_weights = np.hstack([min_weights, max_weights])
    w_LHS = systemLHS * np.sqrt(flat_weights[:, np.newaxis])
    w_RHS = systemRHS * np.sqrt(flat_weights)
    w_lstsq_result = np.linalg.lstsq(w_LHS, w_RHS, rcond=None)[0]

    skew_threshold = 0.20
    is_skewed = delta_median > skew_threshold

    # save everything into a pile of statistics
    xci_stats = {
        'delta': {
            'mean' : float(delta_mean),
            'weighted_mean' : float(weighted_mean),
            'median' : float(delta_median)
        },
        'cluster_center' : {
            'mean' : list(cluster_means),
            'weighted_mean' : list(cluster_weighted_means),
            'median' : list(cluster_medians),
        },
        'skew_threshold' : skew_threshold,
        'linalg_result' : {
            'ratio1' : float(linalg_result[0]),
            'ratio2' : float(1-linalg_result[0]),
            'error' : float(linalg_result[1])
        },
        'lstsq_result' : {
            'ratio1' : float(lstsq_result[0]),
            'ratio2' : float(1-lstsq_result[0]),
            'error' : float(lstsq_result[1])
        },
        'w_lstsq_result' : {
            'ratio1' : float(w_lstsq_result[0]),
            'ratio2' : float(1-w_lstsq_result[0]),
            'error' : float(w_lstsq_result[1])
        },
        'is_skewed' : bool(is_skewed)
    }

    # also save the raw values for visualizations
    data_arrays = {
        'h1_vals' : h1_vals,
        'h2_vals' : h2_vals,
        'h1_unmethyl_reads' : h1_unmethyl_reads,
        'h1_methyl_reads' : h1_methyl_reads,
        'h2_unmethyl_reads' : h2_unmethyl_reads,
        'h2_methyl_reads' : h2_methyl_reads
    }

    return xci_stats, data_arrays

def generateXciImages(xci_stats, xci_data, output_prefix, sample_name):
    '''
    This will generate images that may be useful for XCI visual interpretation
    :param xci_stats: - a dictionary of XCI-related statistics from function `calculateStatistics`
    :param xci_data: - a dictionary of XCI-related numpy arrays from function `calculateStatistics`
    :param output_prefix: - output prefix for saving images to
    :param sample_name: - the provided sample name, mainly for labeling figures
    :return: - None, but will create image files from the `output_prefix`
    '''
    # import for imaging and clustering
    import matplotlib
    matplotlib.use('AGG')
    import matplotlib.pyplot as plt
    import seaborn as sns

    # gather h1 and h2 into some arrays
    h1_vals = xci_data['h1_vals']
    h2_vals = xci_data['h2_vals']
    h1_unmethyl_reads = xci_data['h1_unmethyl_reads']
    h1_methyl_reads = xci_data['h1_methyl_reads']
    h2_unmethyl_reads = xci_data['h2_unmethyl_reads']
    h2_methyl_reads = xci_data['h2_methyl_reads']
    combined_unmethyl_reads = h1_unmethyl_reads + h2_unmethyl_reads
    combined_methyl_reads = h1_methyl_reads + h2_methyl_reads
    
    # statistical calculations
    delta_mean = xci_stats['delta']['mean']
    delta_median = xci_stats['delta']['median']
    skew_threshold = xci_stats['skew_threshold']
    is_skewed = xci_stats['is_skewed']
    ratio1 = xci_stats['linalg_result']['ratio1']
    ratio2 = xci_stats['linalg_result']['ratio2']

    # figure out which of the h1/h2 values is higher
    hap_stack = np.vstack([h1_vals, h2_vals])    
    actual_points = hap_stack.transpose()
    cluster_medians = xci_stats['cluster_center']['median']
    
    # make the .hapscat.png plot
    plt.figure()
    
    # scatter plot the raw data first
    plt.scatter(
        actual_points[:, 0], actual_points[:, 1], # label=cluster_label,
        alpha=0.2, s=5
    )

    # overlay the center-points
    plt.scatter(cluster_medians, cluster_medians[::-1], marker='x', color='k', s=140, label='Median observation')
    plt.scatter([ratio1, ratio2], [ratio2, ratio1], marker='+', color='k', s=140, label='Least-squares ratios')
    
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.plot([0, 1], [0, 1], linestyle='--', color='black')
    plt.xlabel('h1 prob')
    plt.ylabel('h2 prob')
    plt.legend()
    plt.grid()
    plt.title(f'{sample_name}, haplotype methylation ratios\nskew_delta = {delta_mean:.2f}, {delta_median:.2f}; is_skewed({skew_threshold}) = {is_skewed}')
    plt.savefig(f'{output_prefix}.hapscat.png', bbox_inches='tight')
    plt.close()

    # reformat read counts into a heat map
    # original just used the max of all
    # max_value = max(max(unmethyl_reads), max(methyl_reads))

    # take the max of x/y-axis, then average those, ceil it, and then double it
    # impact is that the high points should be center on axes
    max_value = 2 * int(np.ceil(
        np.mean(
            np.max(
                np.vstack([combined_unmethyl_reads, combined_methyl_reads]), 
                axis=0
            )
        )
    ))

    heat_blocks = np.zeros(dtype='int', shape=(max_value+1, max_value+1))
    for i in range(0, len(h1_unmethyl_reads)):
        if h1_unmethyl_reads[i] <= max_value and h1_methyl_reads[i] <= max_value:
            heat_blocks[h1_unmethyl_reads[i], h1_methyl_reads[i]] += 1

        if h2_unmethyl_reads[i] <= max_value and h2_methyl_reads[i] <= max_value:
            heat_blocks[h2_unmethyl_reads[i], h2_methyl_reads[i]] += 1
    
    # make the .read_counts.png plot
    plt.figure()
    ax = sns.heatmap(heat_blocks)
    ax.invert_yaxis()
    plt.xlabel('Unmethylated read count')
    plt.ylabel('Methylated read count')
    plt.title(f'{sample_name}, haplotype-specific read counts')
    plt.savefig(f'{output_prefix}.read_counts.png', bbox_inches='tight')
    plt.close()

def loadPhaseBlocks(phase_block_fn):
    '''
    This will load all chrX phase blocks from a given blocks.tsv file.
    It attempts to auto-detect the phase block source (WhatsHap or HiPhase) while loading.
    :param phase_block_fn: - the phase block .tsv file, can be from either WhatsHap or HiPhase
    :return: - a sorted list of the the phase blocks on chrX as tuples (start, stop, phase_block_id)
    '''
    blocks = []
    fp = open(phase_block_fn, 'r')
    tsv_reader = csv.DictReader(fp, delimiter='\t')
    WHATSHAP_FIELDNAMES = ['#sample', 'chromosome', 'phase_set', 'from', 'to', 'variants']
    if tsv_reader.fieldnames == WHATSHAP_FIELDNAMES:
        phase_source = 'whatshap'
    else:
        phase_source = 'hiphase'

    for row in tsv_reader:
        if phase_source == 'whatshap':
            chrom = row['chromosome']
            start = int(row['from'])
            stop = int(row['to'])
            label = row['phase_set']
        elif phase_source == 'hiphase':
            chrom = row['chrom']
            start = int(row['start'])
            stop = int(row['end'])
            label = row['phase_block_id']
        else:
            raise Exception('should not happen')

        if chrom == 'chrX':
            blocks.append((start, stop+1, label))

    fp.close()
    
    # should generally be sorted already, but lets make it happen
    blocks.sort()
    return blocks

def collapseCpgs(sample_cpgs, phase_block_fn):
    '''
    This will take a selection of CpGs and collapse them into representative phase block CpG.
    Counts get added together as if it was one locus and then methyl probs calculated at the end.
    :param sample_cpgs: - a dictionary of CpG values on a single 
    :param phase_block_fn: - the phase block file to load chrX blocks from
    :return: - a dictionary of CpG values where the keys are now phase block IDs instead of individual CpG loci
    '''
    block_cpgs = {}

    # list of (start, stop, label) tuples, all on chrX
    # TODO: we're hacking to work with just WhatsHap, need to fix this long-term
    blocks = loadPhaseBlocks(phase_block_fn)

    for start_coordinate in sample_cpgs:
        cpg_dict = sample_cpgs[start_coordinate]
        if ('hap1' not in cpg_dict) or ('hap2' not in cpg_dict):
            continue

        block_id = None
        for (start, stop, label) in blocks:
            if start_coordinate >= start and start_coordinate < stop:
                block_id = label
                break

        if block_id == None:
            # it's not _inside_ a phase block, but near one; for now, lets just leave it with the start coordinate
            # ideal solution would have the phase blocks from the source
            block_id = start_coordinate

        # we found it
        if block_id not in block_cpgs:
            block_cpgs[block_id] = {}
        
        # go through each of the methylation categories and add the values
        for k in cpg_dict:
            if k not in block_cpgs[block_id]:
                # this is the first time we hit it, so set to zero initially
                block_cpgs[block_id][k] = {
                    'read_count' : 0,
                    'methyl_count' : 0,
                    'unmethyl_count' : 0,
                    #'methyl_prob': result_prob 
                }
            
            # add them all
            for sub_key in ['read_count', 'methyl_count', 'unmethyl_count']:
                block_cpgs[block_id][k][sub_key] += cpg_dict[k][sub_key]
    
    # add methyl probs for everything
    for block_id in block_cpgs:
        for k in block_cpgs[block_id]:
            block_cpgs[block_id][k]['methyl_prob'] = block_cpgs[block_id][k]['methyl_count'] / block_cpgs[block_id][k]['read_count']

    return block_cpgs

if __name__ == '__main__':
    # set up our parsers
    description = 'Runs a sample through the XCI analysis and reports statistics on XCI status'
    p = ap.ArgumentParser(description=description, formatter_class=ap.RawTextHelpFormatter)

    # optional arguments
    p.add_argument('-i', '--images', dest='images_enabled', action='store_true', default=False, help='enable visualization generation (default: False)')
    p.add_argument('--sample-name', dest='sample_name', default='unspecified', help='set sample name for visualization (default: "unspecified")')
    p.add_argument('--phase-blocks', dest='phase_blocks', default=None, help='phase blocks for grouping observations (default: None)')

    # required args
    p.add_argument('sample_prefix', help='sample prefix from `aligned_bam_to_cpg_scores`')
    p.add_argument('XCI_bed', help='bed file with XCI coordinates')
    p.add_argument('output_prefix', help='the output prefix for all files')
    
    args = p.parse_args()

    # parse the input files
    sample_cpgs = loadSampleXciCpgs(args.sample_prefix, args.XCI_bed)

    # gather statistics and save the results
    xci_statistics, data_arrays = calculateStatistics(sample_cpgs)
    fp = open(f'{args.output_prefix}.xci_summary.json', 'w+')
    json.dump(xci_statistics, fp, indent=4, sort_keys=True)
    fp.close()
    
    # generate images if asked to do so
    if args.images_enabled:
        generateXciImages(xci_statistics, data_arrays, args.output_prefix, args.sample_name)

    if args.phase_blocks != None:
        print('Collapsing results into phase blocks...')
        # phase blocks are provided, also calculate stats at a block level
        block_cpgs = collapseCpgs(sample_cpgs, args.phase_blocks)
        print(f'\tConverted {len(sample_cpgs)} -> {len(block_cpgs)}')

        # save in a separate file
        block_statistics, block_arrays = calculateStatistics(block_cpgs)
        fp = open(f'{args.output_prefix}.phased_xci_summary.json', 'w+')
        json.dump(block_statistics, fp, indent=4, sort_keys=True)
        fp.close()

        #TODO: are images relevant here? my gut says not really
    