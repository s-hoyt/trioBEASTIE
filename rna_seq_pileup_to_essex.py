# Copyright (C)2026 Stephanie H. Hoyt (stephanie.hoyt@duke.edu)
# This is OPEN SOURCE software, released under the MIT License

import sys
import pandas as pd 

het_list = ["0|1", "1|0", "0/1", "1/0"]

#ensembl gtf downloaded from https://jun2026.archive.ensembl.org/info/data/ftp/index.html
def load_ensembl(data_path, version):
    ensembl_gtf = pd.read_csv(data_path + "Homo_sapiens.GRCh38." + version + ".chr.gtf.gz", 
                            sep = "\t", header = None, skiprows = 5,
                            names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
    #remove X, Y, MT annotations
    ensembl_gtf = ensembl_gtf.loc[~ensembl_gtf["seqname"].isin(["X", "Y", "MT"])]
    #make sure remaining chromosome annotations are labeled as ints
    ensembl_gtf["seqname"] = [int(x) for x in ensembl_gtf["seqname"]]
    #just use exons from canonical transcripts
    ensembl_transcripts = ensembl_gtf[(ensembl_gtf["feature"] == "transcript") & ensembl_gtf["seqname"].isin(range(1,23))].reset_index(drop = True)
    is_canon = ["Ensembl_canonical" in x for x in ensembl_transcripts["attribute"]]
    ensembl_transcripts = ensembl_transcripts[is_canon].reset_index(drop=True)
    ensembl_transcripts["transcript_id"] = [x.split(";")[2].split("\"")[1] for x in ensembl_transcripts["attribute"]]

    #want to only use variants in exons
    #need to get all gene ids from ensembl_genes (in attribute column) and then get all exons associated with each gene from ensembl gtf
    ensembl_exons = ensembl_gtf[(ensembl_gtf["feature"] == "exon") & ensembl_gtf["seqname"].isin(range(1,23))].reset_index(drop = True)
    ensembl_exons["transcript_id"] = [x.split(";")[2].split("\"")[1] for x in ensembl_exons["attribute"]]

    return(ensembl_transcripts, ensembl_exons)

#phaser can not handle 'invalid' genotypes, which would involve either
#   a genotyping error
#   or a denovo variant
#so these are removed from the vcf file
def filter_vcf_valid(vcf_path, family_id, indiv_ids):
    child_id, mother_id, father_id = indiv_ids
    trio_vcf = pd.read_csv(vcf_path + family_id + ".wasp.biallelic.het.wheader.vcf", 
                            #compression = 'gzip', 
                            sep = "\t", comment = '#', header = None,
                            names = ["chr", "position", "name", "ref", "alt", 
                                    "qual", "filter", "info", "format", child_id, mother_id, father_id],
                            usecols = ["chr", "position", "ref", "alt", child_id, mother_id, father_id])
    #keep track of which rows will need to be removed
    indices_to_drop = []
    #if at least one parent is homozygous alternate, then child can not be homozygous reference
    indices_to_drop += (list(trio_vcf[((trio_vcf[mother_id] == "1|1") | (trio_vcf[father_id] == "1|1")) & \
                                      (trio_vcf[child_id] == "0|0")].index))
    #if both parents are homozygous ref, then child can not be heterozygous
    indices_to_drop += (list(trio_vcf[((trio_vcf[mother_id] == "0|0") & (trio_vcf[father_id] == "0|0")) & \
                                      ((trio_vcf[child_id] == "0|1") | (trio_vcf[child_id] == "1|0"))].index))
    #if at least one parent is homozygous reference, child can not be homozygous alternate
    indices_to_drop += (list(trio_vcf[((trio_vcf[mother_id] == "0|0") | (trio_vcf[father_id] == "0|0")) & \
                                      (trio_vcf[child_id] == "1|1")].index))
    #if both parents are homozygous alternate, then child can not be heterozygous
    indices_to_drop += (list(trio_vcf[((trio_vcf[mother_id] == "1|1") & (trio_vcf[father_id] == "1|1")) & \
                                      ((trio_vcf[child_id] == "0|1") | (trio_vcf[child_id] == "1|0"))].index))
    return(trio_vcf.drop(indices_to_drop))

def temp_filter_tsv_valid(tsv_path, family_id, indiv_ids):
    homo_ref = ["0|0", "0/0"]
    homo_alt = ["1|1", "1/1"]

    child_id, mother_id, father_id = indiv_ids
    trio_vcf = pd.read_csv(tsv_path + family_id + ".tsv", sep = "\t")
    #keep track of which rows will need to be removed
    indices_to_drop = []
    #if at least one parent is homozygous alternate, then child can not be homozygous reference
    indices_to_drop += (list(trio_vcf[((trio_vcf[mother_id].isin(homo_alt)) | (trio_vcf[father_id].isin(homo_alt))) & \
                                      (trio_vcf[child_id].isin(homo_ref))].index))
    #if both parents are homozygous ref, then child can not be heterozygous
    indices_to_drop += (list(trio_vcf[((trio_vcf[mother_id].isin(homo_ref)) & (trio_vcf[father_id].isin(homo_ref))) & \
                                      (trio_vcf[child_id].isin(het_list))].index))
    #if at least one parent is homozygous reference, child can not be homozygous alternate
    indices_to_drop += (list(trio_vcf[((trio_vcf[mother_id].isin(homo_ref)) | (trio_vcf[father_id].isin(homo_ref))) & \
                                      (trio_vcf[child_id].isin(homo_alt))].index))
    #if both parents are homozygous alternate, then child can not be heterozygous
    indices_to_drop += (list(trio_vcf[((trio_vcf[mother_id].isin(homo_alt)) & (trio_vcf[father_id].isin(homo_alt))) & \
                                      (trio_vcf[child_id].isin(het_list))].index))
    return(trio_vcf.drop(indices_to_drop))


#get the pileups from all the individuals in the trio
#the pileups were generated from the same vcf/positions file, 
#   but different samples have different reads so some positions will be 0 for some individuals even if non zero for others.
#   Therefore, after getting each individual pileup, merge them together, filling missing data with 0 read counts
def merge_pileups(pileup_path, indiv_ids):
    child_id, mother_id, father_id = indiv_ids
    pileups = {}
    for id,indiv in zip([child_id, mother_id, father_id],["child", "father", "mother"]):
        pileup = pd.read_csv(pileup_path + id + ".pileup", 
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
    return pileups_merged

#for a set of SNPs with reads in the trio, want to trim down to sites that are not within a read length of each other
#   so that we are not double counting any reads
#by default do remove triple hets, but add in a parameter to keep them (in case of new model that can handle unphased sites)
def trim(result, indiv_ids, read_len = 200, remove_tri_het = True):
    child_id, mother_id, father_id = indiv_ids
    #label and remove triple hets
    if remove_tri_het:
        tri_het_status = [(result.loc[h, mother_id] in het_list) and \
                            (result.loc[h, father_id] in het_list) and \
                            (result.loc[h, child_id] in het_list) for h in range(len(result))]
        result["trihet"] = tri_het_status
        result = result.loc[~result["trihet"]].reset_index(drop = True)
    if len(result) < 2: #0 or 1 - 0 means it was all triple hets, 1 means we have a phasable site that doesn't need to be trimmed
        return result
    #first find site with deepest total read depth at any heterozygous individual
    deepest_site = 0
    max_depth = 0
    for site in range(len(result)):
        mother = 0
        father = 0 
        child = 0
        if result.loc[site, mother_id] in het_list:
            mother = sum(result.loc[site, ["ref_count_mother", "alt_count_mother"]])
        if result.loc[site, father_id] in het_list:
            father = sum(result.loc[site, ["ref_count_father", "alt_count_father"]])
        if result.loc[site, child_id] in het_list:
            child = sum(result.loc[site, ["ref_count_child", "alt_count_child"]])
        if (mother > max_depth) | (father > max_depth) | (child > max_depth):
            max_depth = max(mother, father, child)
            deepest_site = site
    #then use that as the 'anchor' site for trimming: always keep the anchor site and remove anything within 200 bp of it
    anchor_position = result.loc[deepest_site, "position"]
    #these are all the sites that are NOT the anchor pos, that might need to be trimmed
    result_other = result.loc[result["position"] != anchor_position]
    #these are too close to the anchor pos and do need to be trimmed
    remove = result_other.loc[((result_other["position"] >= anchor_position - read_len) & (result_other["position"] < anchor_position)) \
                                | ((result_other["position"] <= anchor_position + read_len) & (result_other["position"] > anchor_position))]
    # if all sites were too close to the anchor, just return the anchor
    if len(result_other) - len(remove) == 0: 
        return result.loc[[deepest_site,]].reset_index(drop=True)
    else: #if there are some sites that were far enough from the anchor
        #trim them, if there's more than 1 in a given direction
        upstream = result_other.loc[result_other["position"] < anchor_position - read_len].reset_index(drop=True)
        if len(upstream) > 1:
            upstream = trim(upstream, indiv_ids)
        downstream = result_other.loc[result_other["position"] > anchor_position + read_len].reset_index(drop=True)
        if len(downstream) > 1:
            downstream = trim(downstream, indiv_ids)
        #then recombine both directions with anchor
        #this works regardless of which direction (if any) was greater than 1 or if one of them is 0
        return pd.concat([upstream, result.loc[[deepest_site,]], downstream], ignore_index = True)

def make_essex_from_pileup(trio_vcf, pileups_merged, essex_out_dir, family_id, indiv_ids,
                ensembl_transcripts, ensembl_exons, remove = True):
    child_id, mother_id, father_id = indiv_ids
    #GOAL: make essex files with 1) genotypes and 2) pileup allele counts per variant per gene.
    #   allow multiple variants per gene
    
    #make the essex file gene by gene
    file = open(essex_out_dir + family_id + ".essex", "w")

    #all_sites_table = pd.DataFrame()
    for gene in range(len(ensembl_transcripts)):
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
            res_ex = pileup_res_ex.merge(vcf_res_ex, on = ["chr", "position"]) #version for essex
            #res_ex = pileup_res_ex.merge(vcf_res_ex, on = ["chr", "position"], how = "outer") #version for all_sites_table
            res = pd.concat([res, res_ex])
        if len(res) > 0:
            #deduplicate variants in case exon annotations overlap
            res = res.drop_duplicates().reset_index(drop = True)
            res["transcript_id"] = transcript_id
            #save to all_sites table
            ##all_sites_table = pd.concat([all_sites_table, res], ignore_index = True)
            #trim sites to avoid double counting reads
            #this step also checks for triple hets and removes them, which might get rid of a gene if it's all triple hets so need to check len again
            res_trimmed = trim(res, indiv_ids, remove_tri_het = remove)
            if len(res_trimmed) == 0: continue
            file.write("(gene " + transcript_id + "\n")
            for site in range(len(res_trimmed)):
                site_num = res_trimmed.loc[site, "position"]
                #Instead of just putting het variants in the order they are given in vcf (ie either 01 or 10),
                #   always write hets as 01 (as if unphased) because we are going to use the trio phaser 
                #   so vcf phase doesn't matter, and phaser expects 'normalized' variants (no 1 before 0)
                #check if individuals are het
                mother_het = res_trimmed.loc[site, mother_id] in het_list
                father_het = res_trimmed.loc[site, father_id] in het_list
                child_het = res_trimmed.loc[site, child_id] in het_list
                #if individual is heterozygous, set genotype to "0|1" to ensure normalized, else use exact genotype from res
                mother_geno = "0|1" if mother_het else res_trimmed.loc[site, mother_id] 
                father_geno = "0|1" if father_het else res_trimmed.loc[site, father_id] 
                child_geno = "0|1" if child_het else res_trimmed.loc[site, child_id] 
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
    #all_sites_table.to_csv("/hpc/group/majoroslab/stephanie/ceph1463/essex_files_denovo/" + trio + ".rnaseq.tsv", sep = "\t", index = False)

def make_tsv_from_pileup(trio_vcf, pileups_merged, tsv_out_dir, family_id,
                ensembl_transcripts, ensembl_exons):
    #GOAL: tsv files with info about all overlapping variants assigned to their transcripts
    #   even if triple het and even if no overlapping reads, or trimmed for being too close, 
    #   so we have a record of everything
    
    all_sites_table = pd.DataFrame()
    
    for gene in range(len(ensembl_transcripts)):
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
            res_ex = pileup_res_ex.merge(vcf_res_ex, on = ["chr", "position"], how = "outer") #version for all_sites_table
            res = pd.concat([res, res_ex])
        if len(res) > 0:
            #deduplicate variants in case exon annotations overlap
            res = res.drop_duplicates().reset_index(drop = True)
            res["transcript_id"] = transcript_id
            #save to all_sites table
            all_sites_table = pd.concat([all_sites_table, res], ignore_index = True)
    all_sites_table.to_csv(tsv_out_dir + family_id + ".tsv", sep = "\t", index = False)

def make_essex_from_tsv(tsv_out_dir, essex_out_dir, family_id, indiv_ids, remove = True):
    child_id, mother_id, father_id = indiv_ids

    file = open(essex_out_dir + family_id + ".essex", "w")
    sites_df = pd.read_csv(tsv_out_dir + family_id + ".tsv", sep = "\t").dropna().reset_index(drop = True)

    for transcript_id in list(sites_df["transcript_id"].unique()):
        this_transcript = sites_df.loc[sites_df["transcript_id"] ==  transcript_id].dropna().reset_index(drop = True)
        this_transcript = trim(this_transcript, indiv_ids, remove_tri_het = remove)
        if len(this_transcript) == 0: continue
        file.write("(gene " + transcript_id + "\n")
        for site in range(len(this_transcript)):
            site_num = this_transcript.loc[site, "position"]
            mother_het = this_transcript.loc[site, mother_id] in het_list
            father_het = this_transcript.loc[site, father_id] in het_list
            child_het = this_transcript.loc[site, child_id] in het_list
            mother_geno = "0|1" if mother_het else this_transcript.loc[site, mother_id][0] + "|" + this_transcript.loc[site, mother_id][2] 
            father_geno = "0|1" if father_het else this_transcript.loc[site, father_id][0] + "|" + this_transcript.loc[site, father_id][2] 
            child_geno = "0|1" if child_het else this_transcript.loc[site, child_id][0] + "|" + this_transcript.loc[site, child_id][2]
            file.write("\t" + "(site " + str(site_num) + \
                    "  (genotypes (mother " + mother_geno.split("|")[0] + " " + mother_geno.split("|")[1] + \
                            ") (father " + father_geno.split("|")[0] + " " + father_geno.split("|")[1] + \
                                ") (child " + child_geno.split("|")[0] + " " + child_geno.split("|")[1] + ")) \n")
            file.write("\t\t" + \
                " (counts (mother " + str(int(this_transcript.loc[site, "ref_count_mother"])) + " " + str(int(this_transcript.loc[site, "alt_count_mother"])) + \
                    ") (father " + str(int(this_transcript.loc[site, "ref_count_father"])) + " " + str(int(this_transcript.loc[site, "alt_count_father"])) + \
                    ") (child " + str(int(this_transcript.loc[site, "ref_count_child"])) + " " + str(int(this_transcript.loc[site, "alt_count_child"])) + "))) \n")
        file.write(") \n")
        file.flush() 
    file.close()

def main():
    #args = sys.argv[1:]
    # if len(args) != ????:
    #     print("expected arguments are: \n \
    #             <ensembl gtf file> \
    #             <child id> <mother id> <father id>")

    #inputs for load_ensembl
    data_path = "/data/allenlab/stephanie/hg38_reference_files/"
    version = "116"

    #inputs for filter_vcf_valid
    vcf_path = "/data/allenlab/stephanie/sickkids/gs_vcfs/"
    family_id = "CMC_03"
    geno_indiv_ids = [family_id + "_0" + str(i) for i in range(1, 4)]

    #inputs for merge_pileups
    pileup_path = "/data/allenlab/stephanie/sickkids/test/"
    #in order of child, mother, father
    #rna_indiv_ids = ["21_11709", "21_11708", "21_11707"]

    #inputs for make_tsv_from_pileup
    tsv_out_dir = "/data/allenlab/stephanie/sickkids/test/"

    #inputs for make_essex_from_tsv
    essex_out_dir = "/data/allenlab/stephanie/sickkids/test/"

    #call funcs here
    #ensembl_transcripts, ensembl_exons = load_ensembl(data_path, version)
    #filtered_vcf = filter_vcf_valid(vcf_path, family_id, geno_indiv_ids)
    #pileups_merged = merge_pileups(pileup_path, rna_indiv_ids)
    #make_tsv_from_pileup(filtered_vcf, pileups_merged, tsv_out_dir, family_id, ensembl_transcripts, ensembl_exons)
    filtered_tsv = temp_filter_tsv_valid(tsv_out_dir, family_id, geno_indiv_ids)
    filtered_tsv.to_csv(tsv_out_dir + family_id + ".tsv", sep = "\t", index = False)
    make_essex_from_tsv(tsv_out_dir, essex_out_dir, family_id, geno_indiv_ids, remove = True)


if __name__ == "__main__":
    main()
