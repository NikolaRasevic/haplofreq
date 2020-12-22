def haplofreq(peddf):
    """
    
    haplofreq: Software for Determining all Possible Haplotypes and their Estimated Frequencies. 
    
    Description 
    
    haplofreq() is a function that outputs all possible haplotypes as well as their respective frequencies
    for a given inputted genotype dataset. This is accomplished using an EM algorithm developed by Excoffier and 
    Slatkin (1995). 
    
    Usage 
    
    haplofreq(peddf) 
    
    Arguments 
    
    peddf	A .ped file that is read as a Pandas dataframe. The data must not contain any missing genotype data and 
                must only contain bi-allelic SNPs. 
    
    Details 
    
    The algorithm was developed by Excoffier and Slatkin (1995), and takes advantage of the Expectation-Maximization
    algorithm. The expectation step produces a vector of genotype probabilities for all possible genotypes and the
    maximization step consists of estimating the haplotype frequencies for all possible haplotypes using a 
    gene-counting method. The haplofreq function is computationally intensive and it is not recommended that datasets
    containing more than roughly 10 SNPs be inputted into the function.
    
    Output 
    
    A table of all possible haplotypes and haplotype frequencies is outputted. The first column represents the haplotypes
    in square brackets. The second column represents the respective haplotype frequencies estimated by the EM algorithm. 
    
    Example 
    
    data = pd.read_csv("haplofreq\\2snps.ped", delim_whitespace=True, header = None)
    haplofreq(data) 
    """
    import pandas as pd
    import numpy as np
    import itertools
    
    #######
    data = peddf
    data = np.array(data.iloc[:,6:])

    #Reconstructing the ped file
    phenotype = np.zeros(int(len(data.T)/2), dtype = object)

    for i in range(int(len(data.T)/2)):
        phenotype[i] = data[:,int(2*i)] + " " + data[:,int(2*i+1)]
        for x in range(len(phenotype[i])):
            phenotype[i][x] = np.array(phenotype[i][x].split())
        phenotype[i] = np.hstack(phenotype[i])

    phenotype1 = pd.DataFrame(np.zeros([len(phenotype),len(phenotype[0])]))
    for i in range(len(phenotype)):
        phenotype1.iloc[i,:] = phenotype[i]

    phenotype = phenotype1

    #ped file has been reconstructed to phenotype
    #phenotype
    ########
    
    ##############

    #Obtaining sj and cj
    phenotype_code = np.arange(0, len(phenotype.columns)/2, 1)#Phenotype code
    cj = np.zeros(int(len(phenotype.columns)/2))#The number of genotypes that each phenotype contains

    #Loop that assigns values to cj
    for i in range(len(cj)):
        sj = np.sum(phenotype.iloc[:,0-2+2*(i+1)] != phenotype.iloc[:,1-2+2*(i+1)])
        if sj > 0:
            cj[i] = 2**(sj-1)
        if sj == 0:
            cj[i] = 1

    #cj
    ##############

    #############
    #Obtaining all possible haplotypes
    snps = np.zeros(shape=(len(phenotype),2), dtype = object)

    for i in range(len(phenotype)):
        snps[i] = np.unique(phenotype.iloc[i])#Find the major and minor allele of each snps and assign to snps

    snps_input = [0]*len(snps)
    for i in range(len(snps)):
        snps_input[i] = list(snps[i])#Formatting snps so that itertools can find all possible combinations of haplotypes    

    haplotypes = np.array(list(itertools.product(*snps_input))).transpose()#Finding all possible combination of haplotypes
    haplotypes = np.unique(haplotypes, axis = 1)#Subsetting to only unqiue haplotypes
    #haplotypes#All the possible haplotypes
    ##############

    #############
    #Assigning Initial Frequencies to all Haplotypes
    haplotypes_freq = np.ones(len(haplotypes[0]))*(1/len(haplotypes[0]))
    #haplotypes_freq 
    #############

    #############
    #Establising an array of genotypes for each phenotype
    genotype_1 = np.zeros(int(len(phenotype.columns)/2), dtype = object)#alleles of parent 1
    genotype_2 = np.zeros(int(len(phenotype.columns)/2), dtype = object)#alleles of parent 2

    #for i in range(len(genotypes)):
    #    genotypes[i] = np.zeros(shape = (int(len(phenotype)), int(cj[i]*2)), dtype = object)

    geno_snps = np.zeros(shape=(len(phenotype),2), dtype = object)

    for i in range(int(len(phenotype.columns)/2)):#looping through phenotypes
        for t in range(len(geno_snps)):#looping through genotypes
            geno_snps[t][0] = phenotype.iloc[t, 2*i]
            geno_snps[t][1] = phenotype.iloc[t, 2*i+1]
        geno_snps_input = [0]*len(geno_snps)
        for h in range(len(geno_snps)):#looping through genotypes
            geno_snps_input[h] = list(geno_snps[h])
        x = np.array(list(itertools.product(*geno_snps_input))).transpose()
        x = np.unique(x, axis = 1)#x contains all the genotypes of a particular phenotype.
    #    print(x)
        if len(x.T) != 1:#I had to put this condition in so that x values of size 1 can be used.
            #First half of x is assigned to genotype_1
            #Second half of x is assigned to genotype_2
            genotype_1[i] = np.unique(x[:,0:(int(len(x[0])/2))], axis = 1)
            genotype_2[i] = np.unique(x[:,(int(len(x[0])/2)):(int(len(x[0])))], axis = 1)
            genotype_2[i] = np.flip(genotype_2[i], axis = 1)
        if len(x.T) == 1:
            genotype_1[i] = np.unique(x, axis = 1)
            genotype_2[i] = np.unique(x, axis = 1)
            genotype_2[i] = np.flip(genotype_2[i], axis = 1)    

    #Printing the Genotypes
    #print(genotype_1)
    #print(genotype_2)
    ###############

    #####################
    #The EM Alogirithm

    #probability of genotype made up of haplotypes k and l
    def p_hap(k, l):
        k_freq = haplotypes_freq[int(np.where((haplotypes.T == k).all(axis=1))[0])]#finding which haplotypes and their freq == k
        l_freq = haplotypes_freq[int(np.where((haplotypes.T == l).all(axis=1))[0])]#finding which haplotypes and their freq == l
        if k == l:
            p = k_freq**2
        if k != l:
            p = 2*k_freq*l_freq
        return(p)

    #This is the denominator of the expectation equation
    def Pj(genotype_1, genotype_2):
        sum_p = np.zeros(len(genotype_1))
        geno_1 = np.concatenate(genotype_1, axis=1).T
        geno_2 = np.concatenate(genotype_2, axis=1).T
        sum_p_all = np.zeros(len(geno_1))

        for i in range(len(sum_p_all)):
            sum_p_all[i] = p_hap(list(geno_1[i]), list(geno_2[i]))
        count = np.zeros(len(cj))
        count[0] = cj[0]
        for i in range(len(cj)-1):
            count[i+1] = int(count[i] + cj[i+1])
            count = count.astype(int)
        sum_p_all = np.split(sum_p_all, list(count))[0:len(cj)]
        for i in range(len(sum_p)):
            sum_p[i] = np.sum(sum_p_all[i])
        return(sum_p)

    converge = 1

    while converge > 0.000001:

        #Expectation step
        geno_p = np.zeros(len(genotype_1), dtype = object)

        for i in range(len(cj)):
            geno_p[i] = np.zeros(int(cj[i]))

        for i in range(int(len(phenotype.columns)/2)):
            for t in range(len(genotype_1[i].T)):
                #This is the expectation equation
                geno_p[i][t] = (1/(len(phenotype_code)))*(p_hap(list(genotype_1[i][:,t]), (list(genotype_2[i][:,t]))))/(Pj(genotype_1, genotype_2)[i])

        #This is the expectation step for all phenotypes
        #In other words, the probability of each haplotype in the genotypes
        #geno_p

        #Maximization Step

        #Defining the indicator variable
        def δ(haplotype, genotype_1, genotype_2):
            if np.array_equal(haplotype, genotype_1) == True and np.array_equal(haplotype, genotype_2) == True:
                d = 2
            if np.array_equal(haplotype, genotype_1) == False and np.array_equal(haplotype, genotype_2) == False:
                d = 0
            if np.array_equal(haplotype, genotype_1) == True and np.array_equal(haplotype, genotype_2) == False:
                d = 1
            if np.array_equal(haplotype, genotype_1) == False and np.array_equal(haplotype, genotype_2) == True:
                d = 1
            return(d)

        delta_hap = np.zeros(len(haplotypes.T), dtype = object)#Setting a null array(s) in order to obtain indicator variables

        #obtaining values for the inidicator variable
        for x in range(len(delta_hap)):
            delta_p = np.zeros(len(genotype_1), dtype = object)
            for i in range(len(cj)):
                delta_p[i] = np.zeros(int(cj[i]))
            #This is with the 1st haplotype only
            for i in range(int(len(phenotype.columns)/2)):
                for t in range(len(genotype_1[i].T)):
                    delta_p[i][t] = δ(haplotypes[:,x], genotype_1[i][:,t], genotype_2[i][:,t])
            delta_hap[x] = delta_p

        #delta_hap

        #Updating Haplotype Frequencies
        #Using the arrays of indicator variables and haplotype probabilities in respective genotypes
        #this is the maximization step
        haplotypes_freq_new = np.zeros(len(haplotypes_freq))

        for i in range(len(haplotypes_freq)):
            haplotypes_freq_new[i] = (np.sum(np.concatenate((delta_hap[i]*geno_p), axis=0)))/2

        converge = np.max(np.absolute(np.subtract(haplotypes_freq, haplotypes_freq_new)))#Determining the max difference in haplotype freq

        haplotypes_freq = haplotypes_freq_new#Setting the haplotype frequencies to the new one

    #print(haplotypes_freq)
    #print(haplotypes.T)

    output  = np.zeros(shape = (len(haplotypes_freq), 2), dtype = object)
    output[:,0] = list(haplotypes.T)
    output[:,1] = haplotypes_freq
    output = pd.DataFrame(output)
    output.columns = ("Haplotypes", "Frequency")
    return(output)
#####################