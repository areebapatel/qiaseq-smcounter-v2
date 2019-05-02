import sys
from collections import defaultdict

def is_close(a, b, rel_tol=1e-03, abs_tol=1e-03):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def store_outlong_info(outlong):
    ''' Read the outlong file and store information from the file in a dictionary
    :param str outlong: the outlong file from smcounterv2
    :returns information from the outlong file stored in a dictionary
    :rtype dict of dicts
    '''
    OUT_LONG = defaultdict(dict)    
    with open(outlong,'r') as IN:
        for line in IN:
            contents = line.strip('\n').split('\t')        
            if line.startswith('CHROM'):
                header = contents

            chrom,pos,ref,alt = contents[0],contents[1],contents[2],contents[3]
            for i in range(4,len(header)):
                OUT_LONG[chrom+'_'+pos+'_'+ref+'_'+alt][header[i]] = contents[i]
                
    return OUT_LONG

def compare_outlong_info(ORIG_OUT_LONG,NEW_OUT_LONG,verbose):
    ''' Compare outlong info between two smcounter runs
    :param dict of dicts ORIG_OUT_LONG: Original smcounter run being compared to
    :param dict of dicts NEW_OUT_LONG: New smcounter run to test
    :param bool verbose: whether to print information
    :returns (num sites compared, num sites with differing alleles called , num sites with differing metrics called)
    :rtype tuple of ints
    '''
    num_sites = 0
    num_variants_diff = 0
    num_out_differences = 0
    
    out_diff_var = set()
    not_in_new = []
    not_in_old = []
    metric_diff = []
    
    for variant in ORIG_OUT_LONG:
        num_sites+=1
        flag = 0
        if variant not in NEW_OUT_LONG:
            num_variants_diff+=1
            not_in_new.append("Not in NewFile : {var} : {info}".format(var = variant,info = ORIG_OUT_LONG[variant]))

    for variant in NEW_OUT_LONG:
        flag = 0
        if variant not in ORIG_OUT_LONG:
            num_variants_diff+=1
            not_in_old.append("Not in OldFile : {var} : {info}".format(var = variant,info = NEW_OUT_LONG[variant]))
        else:
            for info in NEW_OUT_LONG[variant].keys():
                a = NEW_OUT_LONG[variant][info]
                b = ORIG_OUT_LONG[variant][info]
                try:
                    a = float(a)
                    b = float(b)
                    if not is_close(a, b):
                        flag = 1
                        metric_diff.append(
                            "Different information for variant : {var} : {info} OldFile: {val1} NewFile: {val2}".format(
                                var=variant,info=info,val1=ORIG_OUT_LONG[variant][info],val2=NEW_OUT_LONG[variant][info]))
                except ValueError:
                    try:
                        a = int(a)
                        b = int(a)
                        if not is_close(a, b):
                            flag = 1
                            metric_diff.append("Different information for variant : {var} : {info} OldFile: {val1} NewFile: {val2}".format(
                                var=variant,info=info,val1=ORIG_OUT_LONG[variant][info],val2=NEW_OUT_LONG[variant][info]))
                    except ValueError:
                        if NEW_OUT_LONG[variant][info] != ORIG_OUT_LONG[variant][info]:
                            flag = 1
                            metric_diff.append("Different information for variant : {var} : {info} OldFile: {val1} NewFile: {val2}".format(
                                var=variant,info=info,val1=ORIG_OUT_LONG[variant][info],val2=NEW_OUT_LONG[variant][info]))

            if flag==1:
                num_out_differences+=1

    if verbose:
        for info in not_in_new:
            print info

        for info in not_in_old:
            print info
        for info in metric_diff:
            print info

        print "Total number of sites : {}".format(num_sites)
        print "Number of sites with differing metrics : {}".format(num_out_differences)
        print "Number of sites with differing alleles called  : {}".format(num_variants_diff)

    return (num_sites,num_out_differences,num_variants_diff)

def main(orig_out_long,new_out_long,ci_test,verbose):
    '''
    '''
    ORIG_OUT_LONG = store_outlong_info(orig_out_long)
    NEW_OUT_LONG = store_outlong_info(new_out_long)
    
    num_sites,num_out_differences,num_variants_diff = compare_outlong_info(ORIG_OUT_LONG,NEW_OUT_LONG,verbose)
    
    if ci_test == True:  # If being tested on travis; To do: add more comparisons between the files 
        assert num_variants_diff == 0     # For travis tests, make sure different alleles are not called
        assert num_out_differences <= 10  # metrics can be different slightly (floating point)
    

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print "Run as : python compare_outlong.py <original out long> <new out long> <ci_test: True/False> <verbose: True/False>"
        sys.exit(-1)
    main( sys.argv[1],sys.argv[2],bool(sys.argv[3]),bool(sys.argv[4]))
