# Copyright (C)2026 Stephanie H. Hoyt (stephanie.hoyt@duke.edu)
# This is OPEN SOURCE software, released under the MIT License

import pandas as pd 
from my_functions import trim 
from my_functions import filter_vcf_valid

trios = {"paternal": ["NA12877", "NA12889", "NA12890"], "maternal": ["NA12878", "NA12891", "NA12892"]}
het_list = ["0|1", "1|0"]

ensembl_gtf = pd.read_csv("/hpc/group/majoroslab/stephanie/hg38_reference_files/Homo_sapiens.GRCh38.112.chr.gtf.gz", 
                          sep = "\t", header = None, skiprows = 5,
                          names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
#remove X, Y, MT annotations
ensembl_gtf = ensembl_gtf.loc[~ensembl_gtf["seqname"].isin(["X", "Y", "MT"])]
#make sure remaining chromosome annotations are labeled as ints
ensembl_gtf["seqname"] = [int(x) for x in ensembl_gtf["seqname"]]
#ensembl_genes = ensembl_gtf[(ensembl_gtf["feature"] == "gene") & ensembl_gtf["seqname"].isin(range(1,23))].reset_index(drop = True)
#edit february 11 - just use exons from canonical transcripts
ensembl_transcripts = ensembl_gtf[(ensembl_gtf["feature"] == "transcript") & ensembl_gtf["seqname"].isin(range(1,23))].reset_index(drop = True)
is_canon = ["Ensembl_canonical" in x for x in ensembl_transcripts["attribute"]]
ensembl_transcripts = ensembl_transcripts[is_canon].reset_index(drop=True)
#ensembl_transcripts["gene_id"] = [x.split(";")[0].split("\"")[1] for x in ensembl_transcripts["attribute"]]
ensembl_transcripts["transcript_id"] = [x.split(";")[2].split("\"")[1] for x in ensembl_transcripts["attribute"]]

#Edit November 17 - want to only use variants in exons
#   need to get all gene ids from ensembl_genes (in attribute column) and then get all exons associated with each gene from ensembl gtf
ensembl_exons = ensembl_gtf[(ensembl_gtf["feature"] == "exon") & ensembl_gtf["seqname"].isin(range(1,23))].reset_index(drop = True)
ensembl_exons["transcript_id"] = [x.split(";")[2].split("\"")[1] for x in ensembl_exons["attribute"]]

#Edit Dec 27 - noticed that we have different total number of genes between different runs, which seems to be because code is just taking sooo
#   long to run that it's getting killed at different times but often before it finishes
#Really I should write some code to batch it or something but for now I'm just going to manually check the last gene pick up from there
#lastGene = {"maternal": 45543, "paternal": 45384}

#GOAL: make essex files with 1) genotypes and 2) pileup allele counts per variant per gene.
#   allow multiple variants per gene

#for trio in trios:
for trio in ["maternal"]: #run a separate slurm job for each indiv bc they take so long
    #trio_vcf = filter_vcf_valid(trio)
    #update july 31 2026 - added de novo into phaser, so don't need to filter anymore.
    trio_vcf = pd.read_csv("/hpc/group/majoroslab/stephanie/ceph1463/vcf_" + trio + "_trio/" + trio + "_trio_all_chr.snp.het.vcf.gz", 
                            compression = 'gzip', sep = "\t", comment = '#', header = None,
                            names = ["chr", "position", "name", "ref", "alt", 
                                    "qual", "filter", "info", "format", trios[trio][0], trios[trio][1], trios[trio][2]],
                            usecols = ["chr", "position", "ref", "alt", trios[trio][0], trios[trio][1], trios[trio][2]])
    #get the pileups from all the individuals in the trio
    #the pileups were generated from the same vcf/positions file, 
    #   but different samples have different reads so some positions will be 0 for some individuals even if non zero for others.
    #   Therefore, after getting each individual pileup, merge them together, filling missing data with 0 read counts
    pileups = {}
    for i,indiv in zip(range(3), ["child", "father", "mother"]):
        pileup = pd.read_csv("/cwork/sh700/rna_aligned_encode/GM" + trios[trio][i][2:] + "/align/rep1/" + trio + "_trio_het_sites.pileup", 
                                sep = "\t", header = None,
                                names = ["chr", "position", "ref_allele", "read_depth", "alleles", "scores"])
        ref_count = pileup["alleles"].str.count(r"\.|,")
        pileup["ref_count"] = ref_count
        alt_count = pileup["alleles"].str.count("[acgtACGT]")
        pileup["alt_count"] = alt_count
        total = ref_count + alt_count
        pileup["total_count"] = total
        pileup_noN = pileup[pileup["total_count"] != 0].reset_index(drop = True)
        pileups[indiv] = pileup_noN
    keep_cols = ["chr", "position", "ref_count", "alt_count"]
    #merge all together
    pileups_merged = pileups["child"][keep_cols].merge(pileups["father"][keep_cols].merge(pileups["mother"][keep_cols], 
                                                                                            on = ["chr", "position"], 
                                                                                            suffixes = ["_father", "_mother"], 
                                                                                            how = "outer"), 
                                                        on = ["chr", "position"], 
                                                        how = "outer")
    pileups_merged = pileups_merged.rename(columns = {"ref_count": "ref_count_child", "alt_count": "alt_count_child"})
    #set missing data to 0 read counts
    pileups_merged = pileups_merged.fillna(0)
    #make the essex file gene by gene
    file = open("/hpc/group/majoroslab/stephanie/ceph1463/essex_files_denovo/" + trio + ".rnaseq.essex", "w")
    #Edit January 28 - I need to re-run this to change how I label sites (need to actually know which variants they are!!) and to remove triple hets
    #   but don't want/need to change the genes going into the model
    #So I'm going to cheat and just loop over genes I know came out of the last time I ran this, so that it's faster
    # ensembl_genes_positive = pd.read_csv("/hpc/group/majoroslab/stephanie/ceph1463/essex_files_exons/" + trio + "_genes_w_results.bed", #not a bed file, dunno why i saved it like this
    #                                         sep = "\t", index_col = 0)
    all_sites_table = pd.DataFrame()
    for gene in range(len(ensembl_transcripts)):
    #for gene in range(lastGene[trio], len(ensembl_genes)): #goes with December 27 edit
    #for gene in list(ensembl_genes_positive.index): #goes with January 28 edit
    #for gene in [x for x in list(ensembl_genes_positive.index) if x >= lastGene[trio]]: #combines dec and january edits
        chr = "chr" + str(ensembl_transcripts.loc[gene, "seqname"])
        transcript_id = ensembl_transcripts.loc[gene, "transcript_id"]
        #use all exons, even if the annotations overlap or are duplicated. then de-duplicate variants across all annotatations.
        this_gene = ensembl_exons.loc[ensembl_exons["transcript_id"] == transcript_id].reset_index(drop = True)
        res = pd.DataFrame()
        for ex in range(len(this_gene)):
            start_ex = this_gene.loc[ex, "start"]
            end_ex = this_gene.loc[ex, "end"]
            vcf_res_ex = trio_vcf[(trio_vcf["chr"] == chr) & (trio_vcf["position"] >= start_ex) & (trio_vcf["position"] <= end_ex)]
            pileup_res_ex = pileups_merged[(pileups_merged["chr"] == chr) & (pileups_merged["position"] >= start_ex) & (pileups_merged["position"] <= end_ex)]
            #res_ex = pileup_res_ex.merge(vcf_res_ex, on = ["chr", "position"]) #version for essex
            res_ex = pileup_res_ex.merge(vcf_res_ex, on = ["chr", "position"], how = "outer") #version for all site tsv
            res = pd.concat([res, res_ex])
        if len(res) > 0:
            #deduplicate variants in case exon annotations overlap
            res = res.drop_duplicates().reset_index(drop = True)
            res["transcript_id"] = transcript_id
            #save to all_sites table
            all_sites_table = pd.concat([all_sites_table, res], ignore_index = True)
            #trim sites to avoid double counting reads
            #this step also checks for triple hets and removes them, which might get rid of a gene if it's all triple hets so need to check len again
            res_trimmed = trim(res, trios[trio])
            if len(res_trimmed) == 0: continue
            file.write("(gene " + transcript_id + "\n")
            for site in range(len(res_trimmed)):
                site_num = res_trimmed.loc[site, "position"]
                #Update Aug 8 2024 - Instead of just putting het variants in the order they are given in vcf (ie either 01 or 10),
                #   always write hets as 01 (as if unphased) because we are going to use the trio phaser 
                #   so vcf phase doesn't matter, and phaser expects 'normalized' variants (no 1 before 0)
                #check if individuals are het
                #note to self: there was a bug in the old version of this code that always looked at row 0, not all sites
                #that's why the updated essex files have different variants, not just fewer (which is what I was expecting when I subsetted to just exonic vars)
                mother_het = res_trimmed.loc[site, trios[trio][2]] in het_list
                father_het = res_trimmed.loc[site, trios[trio][1]] in het_list
                child_het = res_trimmed.loc[site, trios[trio][0]] in het_list
                #if individual is heterozygous, set genotype to "0|1" to ensure normalized, else use exact genotype from res
                mother_geno = "0|1" if mother_het else res_trimmed.loc[site, trios[trio][2]] 
                father_geno = "0|1" if father_het else res_trimmed.loc[site, trios[trio][1]] 
                child_geno = "0|1" if child_het else res_trimmed.loc[site, trios[trio][0]] 
                #can change from str(site) to str(site_num) to keep track of how close vars are to each other, but I implemented trimming so we can trust they are at least 200 away
                #changing to site_num not to track how close they are, but so when I do downstream analysis, I know which variants I'm using
                file.write("\t" + "(site " + str(site_num) + \
                        "  (genotypes (mother " + mother_geno.split("|")[0] + " " + mother_geno.split("|")[1] + \
                                ") (father " + father_geno.split("|")[0] + " " + father_geno.split("|")[1] + \
                                    ") (child " + child_geno.split("|")[0] + " " + child_geno.split("|")[1] + ")) \n")
                file.write("\t\t" + \
                    " (counts (mother " + str(int(res_trimmed.loc[site, "ref_count_mother"])) + " " + str(int(res_trimmed.loc[site, "alt_count_mother"])) + \
                        ") (father " + str(int(res_trimmed.loc[site, "ref_count_father"])) + " " + str(int(res_trimmed.loc[site, "alt_count_father"])) + \
                        ") (child " + str(int(res_trimmed.loc[site, "ref_count_child"])) + " " + str(int(res_trimmed.loc[site, "alt_count_child"])) + "))) \n")
            file.write(") \n")
            file.flush() #flush the buffer
    file.close()
    all_sites_table.to_csv("/hpc/group/majoroslab/stephanie/ceph1463/essex_files_denovo/" + trio + ".rnaseq.tsv", sep = "\t", index = False)
