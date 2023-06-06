# MERGES UCSC and ENSEMBL HUMAN GENE ANNOTATIONS for HG19
# author: Vineet Sharma
# Dated : 30th September 2019

import sys
from pprint import pprint as pp
import collections

def gene_search(int_val, int_list, direction):

    """
    This function uses binary search to position any splicing exon at the gene start
    :param int_val: The left value on plus strand and right value in minus strand
    :param int_list: The list of the keys in the dict{min/max-gene_exon: gene_id}
    :return: gene_id
    """
    first_element = 0
    last_element = len(int_list)-1


    while (first_element < last_element) and (last_element - first_element >= 2):
        mid = (first_element + last_element)//2

        if int_val >= int_list[mid]:
            first_element = mid
        else:
            last_element = mid

    if direction == "+":
        return int_list[first_element]
    else:
        return int_list[last_element]

def sorting_tuple(coordinate_tuple_list):

    new_list = []

    for indiv_tup in coordinate_tuple_list:
        a1, b1 = indiv_tup
        if int(a1) > int(b1):
            a1, b1 = b1, a1
        new_list.append((a1, b1))

    return new_list


def genes_expanse(gene_dict):

    gene_coordinate_dict = {}

    for gene, element in gene_dict.items():

        first_element = min(gene_dict[gene]["start"])
        second_element = max(gene_dict[gene]["stop"])
        coordinate_tuple = first_element, second_element
        gene_coordinate_dict.update({gene: coordinate_tuple})

    return gene_coordinate_dict

def find_substituted_exons(exon_list, exns):

    for exn in exns:
        pass


def test_inclusion(exon_list, exns):

    counter = 0

    for exn in exns:
        if exn in exon_list:
            counter += 1
    if counter == len(exon_list):
        return True
    elif counter != 0:
        return False
    else:
        find_substituted_exons(exon_list, exns)


def find_transcripts_containing_miso_exons(gene, iforms, exon_list, gene_dict, splice_category):

    """
    This finds whether the different miso isoforms are represented in the ucsc gencode file.
    Most importantly, considerations include, A/B nomenclature as well as keep these nomenclatures to the fact that
    multiple events of nomenclature may come from the same type or a different type.
    :param iforms:
    :param exon_list:
    :param gene_dict_pointer:
    :return:
    """
    #print("gene", gene)
    #print("iforms", iforms)
    #print("exon_list", exon_list)
    #print("gene_dict", gene_dict)
    #print("splice_category", splice_category)

    for ensgene, ensgene_exnlist in gene_dict.items():
        pass
        if test_inclusion(exon_list, ensgene_exnlist):
            pass

def getting_exon_lists(gene_dict):

    for gene, transcripts in gene_dict.items():
        return gene, transcripts


def get_actual_gene(transcript, put_gene = None, ens_put = None, ens_post = None, ens_prior = None, pre_gene = None, post_gene = None):

    transcript_gene, transcript_tuple = getting_exon_lists(transcript)

    for transcript_exons in transcript_tuple:
        if put_gene != None:
            if transcript_exons in ens_put:
                return put_gene
        elif ens_prior != None:
            if transcript_exons in ens_prior:
                return pre_gene
        elif ens_post != None:
            if transcript_exons in ens_post:
                return post_gene
    return False

def check_other_events_in_dict(event_dict, gene_transcripts):
    event_coord_list = []
    for events in event_dict["before"]:
        before_cord1, before_cord2 = events[0], events[1]
        event_coord_list.append(int(before_cord1))
        event_coord_list.append(int(before_cord2))
    for events in event_dict["after"]:
        after_cord1, after_cord2 = events[0], events[1]
        event_coord_list.append(int(after_cord1))
        event_coord_list.append(int(after_cord2))
    for transcripts, exon_list in gene_transcripts.items():
        for exons in exon_list:
            exn_start, exn_stop = exons
            if int(exn_start) in event_coord_list:
                return True
            else:
                if int(exn_stop) in event_coord_list:
                    return True
    return False

def check_presence_of_ri_intron(check_list, query):
    #print(check_list, query)
    "Determines which one of the events in the query is the one harbouring the ri_intron"

    found = 0
    check_list_exn1, check_list_exn2 = check_list[0]
    new_list = [check_list_exn1, check_list_exn2]
    exon1, exon2 = query[0], query[1]
    exon1_cd1, exon1_cd2 = exon1
    exon2_cd1, exon2_cd2 = exon2

    if len(set(list(exon1))) == 1:
        ri = exon2
        non_ri = exon1
        return ri, non_ri
    if len(set(list(exon2))) == 1:
        ri = exon1
        non_ri=exon2
        return  ri, non_ri

    for exons in query:
        exn1, exn2 = exons
        if (exn1 in new_list) or (exn2 in new_list) :
            found += 1
            try:
                non_ri = exons
            except UnboundLocalError:
                non_ri = None
        else:
            try:
                ri = exons
            except UnboundLocalError:
                ri = None
    if found == 2:
        if int(exon1_cd2) > int(exon2_cd2):
            ri = exon1
            non_ri = exon2
        else:
            ri = exon2
            non_ri = exon1
    return ri, non_ri

def make_gene_annotation_for_writing(fh, chr, gn, gd, transcript = None, exon_list = None, begin = None,
                                     end = None, annot = None, feature = None):
    events_list = []
    if annot != None and feature == "gene":
        for transcript, event_dict in annot.items():
            for event, loci in event_dict.items():
                events_list.append(event)
        event_strg = "|".join(list(set(events_list)))
        line = chr+"\t"+"Ucsc_Ensmbl_merged"+"\t"+feature+"\t"+str(begin)+"\t"+str(end)+"\t"+"."+"\t"+gd+"\t"+"."+\
               "\t"+"gene_id "+ "\"" +gn+ "\""+"; "+"events_in_gene "+event_strg+";"+"\n"
        fh.write(line)

    if feature == "transcript" and annot != None:
        #print(transcript)
        transcript_events_list = []
        events = annot[transcript]
        #print(events)
        new_event = str()
        for event, loci_list in events.items():
            try:
                coordinates_present = "#".join(loci_list)
                new_event = new_event+str(event)+": "+coordinates_present+"|"
            except TypeError:
                new_event = new_event+str(event)+": "+str(loci_list)+"|"
        new_line = chr+"\t"+"Ucsc_Ensmbl_merged"+"\t"+feature+"\t"+str(begin)+"\t"+str(end)+"\t"+"."+"\t"\
                   + gd + "\t" + "." + "\t" + "gene_id " + "\"" +gn+"\""+"; "+"transcript_id "+ "\"" + transcript + "\"" + "; "+"events_in_transcripts "+ new_event
        line = new_line[:-1]+";"+"\n"
        fh.write(line)
        for exon in exon_list:
            exon_start = exon[0]
            exon_stop = exon[1]
            line = chr + "\t" + "Ucsc_Ensmbl_merged" + "\t" + "exon" + "\t" + str(exon_start) + "\t" + str(exon_stop) + "\t" + "." + "\t" + gd + "\t" + "." + "\t" + "gene_id " + "\"" + gn + "\"" + "; " + "transcript_id "+ "\"" +transcript+ "\"" +";"+"\n"
            fh.write(line)

    if annot == None and feature == "gene":
        line = chr + "\t" + "Ucsc_Ensmbl_merged" + "\t" + feature + "\t" + str(begin) + "\t" + str(end) + "\t" + "." + "\t" + gd + "\t" + "." + "\t" + "gene_id " + "\"" + gn + "\"" + "; " + "\n"
        fh.write(line)
    if annot == None and feature == "transcript":
        line = chr + "\t" + "Ucsc_Ensmbl_merged" + "\t" + feature + "\t" + str(begin) + "\t" + str(end) + "\t" + "." + "\t" + gd + "\t" + "." + "\t" + "gene_id " + "\"" + gn + "\"" + "; " + "transcript_id " + "\"" +transcript + "\"" + ";" + "\n"
        fh.write(line)
        for exon in exon_list:
            exon_start = exon[0]
            exon_stop = exon[1]
            line = chr + "\t" + "Ucsc_Ensmbl_merged" + "\t" + "exon" + "\t" + str(exon_start) + "\t" + str(exon_stop) + "\t" + "." + "\t" + gd + "\t" + "." + "\t" + "gene_id " + "\"" +gn + "\"" + "; " + "transcript_id " + "\"" + transcript+ "\"" + ";" + "\n"
            fh.write(line)


class ReadGTFFile:

    def __init__(self, gtf_fn):

        self.gtf_fn = gtf_fn
        self.line_list = self.read_gtf()

    def read_gtf(self):

        with open(self.gtf_fn) as file:
            return self.all_rows_list(file)

    def line_generator(self, src):

        for line in src.readlines():
            lines = line.split("\t")
            if len(lines[0].split("_")) == 1:
                yield lines

    def all_rows_list(self, src):

        line_list = []

        for lines in self.line_generator(src):
            line_list.append(lines)

        return line_list


class EnsemblGTFs(ReadGTFFile):

    """
    class to parse EnsemblGTFs
    """

    def __init__(self, En_fn):

        self.En_fn = En_fn
        ReadGTFFile.__init__(self, self.En_fn)
        self.gene_id_transcript_id_exon_tuple_dict, self.tr_dict_EN, self.ens_gene_id_chr_dict, \
        self.gene_id_direction_dict, self.gene_id_exon_tuple_dict = self.make_dict_En()

    def make_dict_En(self):

        """
        This maps transcripts to genes from the Ensembl GTF for eg. 'ENSMUST00000191675': 'ENSMUSG00000102668'
        :return:dict of transcript to genes
        """

        transcript_gene_id_dict = {}
        gene_id_chromosome_dict = {}
        gene_id_transcript_id_exon_tuple_dict = {}
        gene_id_direction_dict = {}
        gene_id_exon_tuple_dict = {}

        for line in self.line_list:
            if line[2] == "exon":
                direction = line[6]
                chromosome = line[0]
                exon_tuple = line[3], line[4]
                for description in line[8].split(";"):
                    if description.lstrip().startswith("gene_id"):
                        gene_id = description.split(" ")[1][1:-1]
                    if description.lstrip().startswith("transcript_id"):
                        transcript_id = description.lstrip().split(" ")[1][1:-1]
                if not gene_id in gene_id_transcript_id_exon_tuple_dict.keys():
                    gene_id_exon_tuple_dict[gene_id] = []
                    gene_id_exon_tuple_dict[gene_id].append(exon_tuple)
                    gene_id_direction_dict[gene_id] = direction
                    gene_id_transcript_id_exon_tuple_dict[gene_id] = {}
                    gene_id_transcript_id_exon_tuple_dict[gene_id][transcript_id] = []
                    gene_id_transcript_id_exon_tuple_dict[gene_id][transcript_id].append(exon_tuple)
                    transcript_gene_id_dict[transcript_id] = gene_id
                    gene_id_chromosome_dict[gene_id] = chromosome
                else:
                    gene_id_exon_tuple_dict[gene_id].append(exon_tuple)
                    if not transcript_id in gene_id_transcript_id_exon_tuple_dict[gene_id].keys():
                        gene_id_transcript_id_exon_tuple_dict[gene_id][transcript_id] = []
                        gene_id_transcript_id_exon_tuple_dict[gene_id][transcript_id].append(exon_tuple)
                        transcript_gene_id_dict[transcript_id] = gene_id
                    else:
                        gene_id_transcript_id_exon_tuple_dict[gene_id][transcript_id].append(exon_tuple)

        return gene_id_transcript_id_exon_tuple_dict, transcript_gene_id_dict, gene_id_chromosome_dict, \
               gene_id_direction_dict, gene_id_exon_tuple_dict


class UCSCGencodeComprehensive(ReadGTFFile):

    """
    class to parse the UCSC Gencode file
    """

    def __init__(self, ug_fn):
        self.ug_fn = ug_fn
        ReadGTFFile.__init__(self, self.ug_fn)
        self.ucsc_genedict, self.ucsc_gene_direction_dict, self.ucsc_gene_id_chr_dict = self.make_dict_ug()
        #pp(self.ucsc_genedict)

    def make_dict_ug(self):

        """
        Makes dicts of the gene:transcript:exon_tuples
        :return:
        """
        gene_direction_dict = {}
        gene_dict = {}
        gene_id_chromosome_dict = {}

        for lines in self.line_list:

            if lines[2] == "exon":

                exon_tuple = lines[3], lines[4]
                for description in lines[8].split(";"):
                    if description.lstrip().startswith("gene_id"):
                        gene_id = description.split(" ")[1].split(".")[0][1:]
                        #gene = gene_id.split(".")
                        #gene_id = (gene[0][1:])
                    if description.lstrip().startswith("transcript_id"):
                        transcript_id = description.split(" ")[2].split(".")[0][1:]
                        #transcript = transcript_id.split(".")
                        #transcript_id = transcript[0][1:]
                if gene_id not in gene_dict.keys():
                    gene_id_chromosome_dict[gene_id] = lines[0]
                    gene_direction_dict[gene_id] = lines[6]
                    gene_dict[gene_id] = {}
                    gene_dict[gene_id][transcript_id] = []
                    gene_dict[gene_id][transcript_id].append(exon_tuple)
                else:
                    if transcript_id not in gene_dict[gene_id].keys():
                        gene_dict[gene_id][transcript_id] = []
                        gene_dict[gene_id][transcript_id].append(exon_tuple)
                    else:
                        gene_dict[gene_id][transcript_id].append(exon_tuple)

        return gene_dict, gene_direction_dict, gene_id_chromosome_dict

class SplicingMISOFeatures(ReadGTFFile):

    """
    Parses MISO annotations
    """

    def __init__(self, miso_fn):
        self.miso_fn = miso_fn
        self.miso_geneid_dict = {}
        ReadGTFFile.__init__(self, self.miso_fn)
        self.miso_genedict, self.miso_gene_direction_dict, self.miso_gene_id_chr_dict, self.gene_coordinate_search = \
            self.make_dict_MI()
        #pp(self.miso_geneid_dict["RI"])

    def make_dict_MI(self):

        gene_dict = {}
        gene_direction_dict = {}
        gene_id_chromosome_dict = {}
        gene_coordinate_look_dict = {}

        for lines in self.line_list:
            exon_tuple = lines[3], lines[4]
            chromosome = lines[0]
            miso_class = lines[1]
            for description in lines[8].split(";"):
                if description.lstrip().startswith("gene_id"):
                    gene_id = description.lstrip().split(" ")[1]
                if description.lstrip().startswith("transcript_id"):
                    transcript_id = str(description.lstrip().split(" ")[1][-2:-1])
            if chromosome not in gene_dict.keys():
                gene_dict[chromosome] = {}
            if gene_id not in gene_dict[chromosome].keys():
                try:
                    if gene_direction == "+":
                        coordinate_to_look = min(exon_stop_list)
                        gene_coordinate_look_dict[old_gene_id] = coordinate_to_look
                    else:
                        coordinate_to_look = max(exon_start_list)
                        gene_coordinate_look_dict[old_gene_id] = coordinate_to_look
                except UnboundLocalError:
                    pass

                gene_direction = lines[6]
                gene_direction_dict[gene_id] = gene_direction
                gene_dict[chromosome][gene_id] = {}
                gene_dict[chromosome][gene_id][transcript_id] = []
                gene_dict[chromosome][gene_id][transcript_id].append(exon_tuple)
                exon_start_list = []
                exon_stop_list = []
                exon_start_list.append(int(lines[3]))
                exon_stop_list.append(int(lines[4]))
                gene_id_chromosome_dict[gene_id] = chromosome
                self.make_miso_coordinate_dict_from_geneid(gene_id, gene_direction, miso_class, chromosome)
            else:
                old_gene_id = gene_id
                exon_start_list.append(int(lines[3]))
                exon_stop_list.append(int(lines[4]))
                if transcript_id not in gene_dict[chromosome][gene_id].keys():
                    gene_dict[chromosome][gene_id][transcript_id] = []
                    gene_dict[chromosome][gene_id][transcript_id].append(exon_tuple)
                else:
                    gene_dict[chromosome][gene_id][transcript_id].append(exon_tuple)
        if gene_direction == "+":
            coordinate_to_look = min(exon_stop_list)
            gene_coordinate_look_dict[old_gene_id] = coordinate_to_look
        else:
            coordinate_to_look = max(exon_start_list)
            gene_coordinate_look_dict[old_gene_id] = coordinate_to_look


        return gene_dict, gene_direction_dict, gene_id_chromosome_dict, gene_coordinate_look_dict

    def make_miso_coordinate_dict_from_geneid(self, gene_id, gene_direction, miso_class, chromosome):

        if miso_class == "A3SS":
            if "A3SS" not in self.miso_geneid_dict.keys():
                self.miso_geneid_dict["A3SS"] = {}
            if gene_id not in self.miso_geneid_dict["A3SS"].keys():
                self.miso_geneid_dict["A3SS"][gene_id] = {}
                self.miso_geneid_dict["A3SS"][gene_id].update({"before":[], "after":[], "important_coordinates":[]})
                                                # before & after are only relative to MISO nomenclature not on sequence
                coordinate_elements = gene_id.split("@")
                for coordinate_structure in coordinate_elements:
                    check_a3ss = coordinate_structure.split(":")
                    if "|" in check_a3ss[1]:
                        f1, s1 = check_a3ss[1].split("|")
                        self.miso_geneid_dict["A3SS"][gene_id]["important_coordinates"].append((f1, s1))
                        self.miso_geneid_dict["A3SS"][gene_id]["after"].append((f1, check_a3ss[2]))
                        self.miso_geneid_dict["A3SS"][gene_id]["after"].append((s1, check_a3ss[2]))
                    else:
                        self.miso_geneid_dict["A3SS"][gene_id]["before"].append((check_a3ss[1], check_a3ss[2]))
        if miso_class == "A5SS":
            if "A5SS" not in self.miso_geneid_dict.keys():
                self.miso_geneid_dict["A5SS"] = {}
            if gene_id not in self.miso_geneid_dict["A5SS"].keys():
                self.miso_geneid_dict["A5SS"][gene_id] = {}
                self.miso_geneid_dict["A5SS"][gene_id].update({"before":[], "after":[], "important_coordinates":[]})
                coordinate_structure = gene_id.split("@")
                for coordinate_elements in coordinate_structure:
                    check_a5ss = coordinate_elements.split(":")
                    if "|" in check_a5ss[2]:
                        f1, s1 = check_a5ss[2].split("|")
                        self.miso_geneid_dict["A5SS"][gene_id]["important_coordinates"].append((f1, s1))
                        self.miso_geneid_dict["A5SS"][gene_id]["before"].append((check_a5ss[1], f1))
                        self.miso_geneid_dict["A5SS"][gene_id]["before"].append((check_a5ss[1], s1))
                    else:
                        self.miso_geneid_dict["A5SS"][gene_id]["after"].append((check_a5ss[1], check_a5ss[2]))

        if miso_class == "MXE":
            if "MXE" not in self.miso_geneid_dict.keys():
                self.miso_geneid_dict["MXE"] = {}
            if gene_id not in self.miso_geneid_dict["MXE"].keys():
                self.miso_geneid_dict["MXE"][gene_id] = {}
                self.miso_geneid_dict["MXE"][gene_id].update({"before":[], "after":[], "important_coordinates":[]})
                coordinate_structure = gene_id.split("@")
                chrmsm1, start1, stop1, direction1 = coordinate_structure[0].split(":")
                chrmsm2, start2, stop2, direction2 = coordinate_structure[1].split(":")
                chrmsm3, start3, stop3, direction3 = coordinate_structure[2].split(":")
                chrmsm4, start4, stop4, direction4 = coordinate_structure[3].split(":")
                self.miso_geneid_dict["MXE"][gene_id]["important_coordinates"].append((start2, stop2))
                self.miso_geneid_dict["MXE"][gene_id]["important_coordinates"].append((start3, stop3))
                self.miso_geneid_dict["MXE"][gene_id]["before"].append((start1, stop1))
                self.miso_geneid_dict["MXE"][gene_id]["after"].append((start4, stop4))


        if miso_class == "RI":
            if "RI" not in self.miso_geneid_dict.keys():
                self.miso_geneid_dict["RI"] = {}
            if gene_id not in self.miso_geneid_dict["RI"].keys():
                coordinates = []
                self.miso_geneid_dict["RI"][gene_id] = {}
                self.miso_geneid_dict["RI"][gene_id].update({"before":[], "after":[], "important_coordinates":[]})
                """
                For RI, there is nothing like before/after but I am distributing the nonRI and the RI version which 
                starts upstream as the before and later one as after
                """
                coordinate_structure = gene_id.split("@")
                for coordinate_elements in coordinate_structure:
                    events = coordinate_elements.split(":")
                    event_list = events[1].split("-")
                    for event in event_list:
                        coordinates.append(int(event))
                coordinates.sort()
                self.miso_geneid_dict["RI"][gene_id]["important_coordinates"].append((coordinates[1], coordinates[2]))
                self.miso_geneid_dict["RI"][gene_id]["before"].append((coordinates[0], coordinates[1]))
                self.miso_geneid_dict["RI"][gene_id]["before"].append((coordinates[0], coordinates[3]))
                self.miso_geneid_dict["RI"][gene_id]["after"].append((coordinates[2],coordinates[3]))


        if miso_class == "SE":
            if "SE" not in self.miso_geneid_dict.keys():
                self.miso_geneid_dict["SE"] = {}
            if gene_id not in self.miso_geneid_dict["SE"].keys():
                self.miso_geneid_dict["SE"][gene_id] = {}
                self.miso_geneid_dict["SE"][gene_id].update({"before":[], "after":[], "important_coordinates":[]})
                coordinate_structure = gene_id.split("@")
                coordinate_element1 = coordinate_structure[0].split(":")
                coordinate_element2 = coordinate_structure[1].split(":")
                coordinate_element3 = coordinate_structure[2].split(":")
                self.miso_geneid_dict["SE"][gene_id]["before"].append((coordinate_element1[1], coordinate_element1[2]))
                self.miso_geneid_dict["SE"][gene_id]["important_coordinates"].append((coordinate_element2[1],
                                                                                      coordinate_element2[2]))
                self.miso_geneid_dict["SE"][gene_id]["after"].append((coordinate_element3[1], coordinate_element3[2]))

    def sorting_miso_geneid_dict(self):

        for splicetype, gene_dict in self.miso_geneid_dict.items():
            for gene, coordinate_dicts in gene_dict.items():
                for relative_coordinate_placement, coordinate_tuple_lists in coordinate_dicts.items():
                    sorted_list = sorting_tuple(coordinate_tuple_lists)
                    self.miso_geneid_dict[splicetype][gene].update({relative_coordinate_placement: sorted_list})


class CheckOverlaps(UCSCGencodeComprehensive, SplicingMISOFeatures, EnsemblGTFs):

    """
    Things to note are: 1. Though if there have been assembling of genes that encompass the UCSC + Ensembl + MISO
    annotation but  in the end it is only the UCSC+Ensembl annotation that has been reported in the final gtf to make
    sure that we are not annotating from our own.
    """

    def __init__(self, ug_fn, miso_fn, En_fn, merged_fn):
        self.testing_counter = 0
        self.ug_fn = ug_fn
        self.miso_fn = miso_fn
        self.En_fn = En_fn
        self.merged_fn = merged_fn 
        EnsemblGTFs.__init__(self, self.En_fn)
        UCSCGencodeComprehensive.__init__(self, self.ug_fn)
        SplicingMISOFeatures.__init__(self, self.miso_fn)
        self.gene_dict_with_MISO_features = {}
        self.ensmug_genes, self.ucsc_genes, self.gene_direction, \
        self.ucsc_genes_converted =  self.generate_new_dicts_for_UCSCGencode()
        self.genes_converted_to_ensemble_genes, self.ucsc_specific_genes, self.gene_start_dict, \
        self.ensembl_gene_coordinate_dict = self.gene_boundaries()
        self.final_dict, self.ucsc_genes_to_be_ignored = self.combine_ucsc_specific_genes_to_ensembl_gene_dicts()
        self.gene_start_gene_id_dict, self.gene_start_chr_dict_list = self.make_new_boundaries()
        self.final_miso_dict = self.tester()
        self.gene_boundaries, self.transcript_boundaries = self.determine_start_stop_gene_boundaries_for_gtf()
        self.final_miso_annotation_dict = self.squeeze_unique_events()

    def generate_new_dicts_for_UCSCGencode(self):

        genes_converted_to_ENSMUG = {}
        genes_unable_to_converted_to_ENSMUG = {}
        genes_direction_dict = {}
        ucsc_gene_converted = []


        # Here I realise that there are coordinates overlaps on the UCSC annotation but not on the ensembl annotation.
        # I chose to take the ensembl annotation as that is the base for the annotation and ignore the UCSC annotation
        # for the transcript in the gene

        for gene, transcript_dict in self.ucsc_genedict.items():
            ucsc_gene = gene
            try:
                ens_gene = self.tr_dict_EN[ucsc_gene]
                if not ens_gene in genes_converted_to_ENSMUG.keys():
                    genes_direction_dict[ens_gene] = self.ucsc_gene_direction_dict[gene]
                    genes_converted_to_ENSMUG[ens_gene] = {}
                    genes_converted_to_ENSMUG[ens_gene].update({gene: self.gene_id_transcript_id_exon_tuple_dict[ens_gene][gene]})
                    ucsc_gene_converted.append(ucsc_gene)
                else:
                    genes_converted_to_ENSMUG[ens_gene].update({gene: self.gene_id_transcript_id_exon_tuple_dict[ens_gene][gene]})
                    ucsc_gene_converted.append(ucsc_gene)
            except KeyError:
                if not gene in genes_unable_to_converted_to_ENSMUG.keys():
                    genes_direction_dict[gene] = self.ucsc_gene_direction_dict[gene]
                    genes_unable_to_converted_to_ENSMUG[gene] = {}
                    genes_unable_to_converted_to_ENSMUG[gene].update(transcript_dict)

        # These genes not present in the UCSC gtf are to be incorporated in the genes_converted_to_ENSMUG dict
        # because these genes still have a ENSMUG nomenclature whereas genes_unable_to_converted_to_ENSMUG
        # will have ENSMUST nomenclature and have to be dealt separately.

        for gene_id, transcript_id_dict in self.gene_id_transcript_id_exon_tuple_dict.items():
            if gene_id not in genes_converted_to_ENSMUG.keys():
                genes_converted_to_ENSMUG[gene_id] = transcript_id_dict
                genes_direction_dict[gene_id] = self.gene_id_direction_dict[gene_id]

        return genes_converted_to_ENSMUG, genes_unable_to_converted_to_ENSMUG, genes_direction_dict, ucsc_gene_converted

    def gene_boundaries(self):

        gene_exon_tuples_dict = {}
        ucsc_transcript_tuples_dict = {}

        for genes, transcript_dicts in self.ensmug_genes.items():
            for transcripts, exon_tuple_list in transcript_dicts.items():
                for exon_tuples in exon_tuple_list:
                    if not genes in gene_exon_tuples_dict.keys():
                        gene_exon_tuples_dict[genes] = {}
                        gene_exon_tuples_dict[genes]["start"] = set()
                        gene_exon_tuples_dict[genes]["stop"] = set()
                        gene_exon_tuples_dict[genes]["start"].add(int(exon_tuples[0]))
                        gene_exon_tuples_dict[genes]["stop"].add(int(exon_tuples[1]))
                    else:
                        gene_exon_tuples_dict[genes]["start"].add(int(exon_tuples[0]))
                        gene_exon_tuples_dict[genes]["stop"].add(int(exon_tuples[1]))

        for genes, transcript_dicts in self.ucsc_genes.items():
            for transcripts, exon_tuple_list in transcript_dicts.items():
                for exon_tuples in exon_tuple_list:
                    if not genes in ucsc_transcript_tuples_dict.keys():
                        ucsc_transcript_tuples_dict[genes] = {}
                        ucsc_transcript_tuples_dict[genes]["start"] = set()
                        ucsc_transcript_tuples_dict[genes]["stop"] = set()
                        ucsc_transcript_tuples_dict[genes]["start"].add(int(exon_tuples[0]))
                        ucsc_transcript_tuples_dict[genes]["stop"].add(int(exon_tuples[1]))
                    else:
                        ucsc_transcript_tuples_dict[genes]["start"].add(int(exon_tuples[0]))
                        ucsc_transcript_tuples_dict[genes]["stop"].add(int(exon_tuples[1]))

        genes_converted_to_ensemble_genes = genes_expanse(gene_exon_tuples_dict)
        ucsc_specific_genes = genes_expanse(ucsc_transcript_tuples_dict)

        ensembl_gene_coordinate_dict = {}
        gene_start_dict = {}
        gene_start_dict_coordinates_to_avoid = {}

        for gene, coordinate_tuple in genes_converted_to_ensemble_genes.items():
            c1, c2 = coordinate_tuple
            chrmsm = self.ens_gene_id_chr_dict[gene]
            if not chrmsm in ensembl_gene_coordinate_dict.keys():
                ensembl_gene_coordinate_dict[chrmsm] = []
            if self.gene_direction[gene] == "+":
                if not chrmsm in gene_start_dict.keys():
                    gene_start_dict[chrmsm] = {}
                    gene_start_dict[chrmsm][int(coordinate_tuple[0])] = gene
                    ensembl_gene_coordinate_dict[chrmsm].append(int(coordinate_tuple[0]))
                else:
                    if int(coordinate_tuple[0]) not in gene_start_dict[chrmsm].keys():
                        gene_start_dict[chrmsm][int(coordinate_tuple[0])] = gene
                    else:
                        old_gene = gene_start_dict[chrmsm][int(coordinate_tuple[0])]
                        if type(old_gene) == list:
                            gene_start_dict[chrm][int(coordinate_tuple[0])].append(gene)
                        else:
                            del gene_start_dict[chrmsm][int(coordinate_tuple[0])]
                            gene_start_dict[chrmsm][int(coordinate_tuple[0])] = [old_gene, gene]
                    if not int(coordinate_tuple[0]) in ensembl_gene_coordinate_dict[chrmsm]:
                        ensembl_gene_coordinate_dict[chrmsm].append(int(coordinate_tuple[0]))
            else:
                if not chrmsm in gene_start_dict.keys():
                    gene_start_dict[chrmsm] = {}
                    gene_start_dict[chrmsm][int(coordinate_tuple[1])] = gene
                    ensembl_gene_coordinate_dict[chrmsm].append(int(coordinate_tuple[1]))
                else:
                    if int(coordinate_tuple[1]) not in gene_start_dict[chrmsm].keys():
                        gene_start_dict[chrmsm][int(coordinate_tuple[1])] = gene
                    else:
                        old_gene = gene_start_dict[chrmsm][int(coordinate_tuple[1])]
                        if type(old_gene) == list:
                            gene_start_dict[chrm][int(coordinate_tuple[1])].append(gene)
                        else:
                            del gene_start_dict[chrmsm][int(coordinate_tuple[1])]
                            gene_start_dict[chrmsm][int(coordinate_tuple[1])] = [old_gene, gene]
                    if not int(coordinate_tuple[0]) in ensembl_gene_coordinate_dict[chrmsm]:
                        ensembl_gene_coordinate_dict[chrmsm].append(int(coordinate_tuple[1]))

        for chrm, coordn_list in ensembl_gene_coordinate_dict.items():
            ensembl_gene_coordinate_dict[chrm].sort()

        return genes_converted_to_ensemble_genes, ucsc_specific_genes, \
               gene_start_dict, ensembl_gene_coordinate_dict

    def combine_ucsc_specific_genes_to_ensembl_gene_dicts(self):

        """
        One of the problems identified is that the start of the transcript to be annotated to the ensembl gene data base
        could have the start earlier than the ensembl gene boundaries. It needs to be checked for both directions.
        :return:
        """
        ucsc_transcripts_to_be_ignored = []
        final_dict = self.ensmug_genes

        i = 0

        for genes, transcripts in self.ucsc_specific_genes.items():
            actual_gene = False

            i += 1
            direction = self.ucsc_gene_direction_dict[genes]
            if direction == "+":
                gene_start = transcripts[0]
            else:
                gene_start = transcripts[1]

            chrmsm = self.ucsc_gene_id_chr_dict[genes]

            gene_putative_location = gene_search(int(gene_start), self.ensembl_gene_coordinate_dict[chrmsm], direction)

            # This searches for the putative gene location in a list of sorted starting coordinates

            gene_putative_index = self.ensembl_gene_coordinate_dict[chrmsm].index(gene_putative_location)

            # This gives back the index of the gene_putative_location. Further the test happens for the transcripts
            # that might be beyond extremities of teh gene coordinates and thus may give a false location and hence
            # false gene nomenclature
            
            if type(self.gene_start_dict[chrmsm][self.ensembl_gene_coordinate_dict[chrmsm][gene_putative_index]]) is not list:
                putative_gene = self.gene_start_dict[chrmsm][self.ensembl_gene_coordinate_dict[chrmsm][gene_putative_index]]

                if (gene_putative_index != 0) and (gene_putative_index != len(self.ensembl_gene_coordinate_dict[chrmsm]) - 1):
                    if type(self.gene_start_dict[chrmsm][self.ensembl_gene_coordinate_dict[chrmsm][gene_putative_index - 1]]) is not list:
                            prior_gene = self.gene_start_dict[chrmsm][self.ensembl_gene_coordinate_dict[chrmsm][gene_putative_index - 1]]
                            ens_gene_prior = self.gene_id_exon_tuple_dict[prior_gene]
                    else:
                        prior_gene = None
                        ens_gene_prior = None

                    if type(self.gene_start_dict[chrmsm][self.ensembl_gene_coordinate_dict[chrmsm][gene_putative_index + 1]]) is not list:
                        posterior_gene = self.gene_start_dict[chrmsm][self.ensembl_gene_coordinate_dict[chrmsm][gene_putative_index + 1]]
                        ens_gene_post = self.gene_id_exon_tuple_dict[posterior_gene]
                    else:
                        posterior_gene = None
                        ens_gene_post = None

                    actual_gene = get_actual_gene(self.ucsc_genes[genes], put_gene = putative_gene,
                                                  ens_put = self.gene_id_exon_tuple_dict[putative_gene],
                                                  ens_prior = ens_gene_prior, ens_post = ens_gene_post,
                                                  pre_gene = prior_gene, post_gene = posterior_gene)
            else:
                putative_gene = None
                if (gene_putative_index != 0) and (gene_putative_index != len(self.ensembl_gene_coordinate_dict[chrmsm]) - 1):
                    if type(self.gene_start_dict[chrmsm][self.ensembl_gene_coordinate_dict[chrmsm][gene_putative_index - 1]]) is not list:
                            prior_gene = self.gene_start_dict[chrmsm][self.ensembl_gene_coordinate_dict[chrmsm][gene_putative_index - 1]]
                            ens_gene_prior = self.gene_id_exon_tuple_dict[prior_gene]
                    else:
                        prior_gene = None
                        ens_gene_prior = None

                    if type(self.gene_start_dict[chrmsm][self.ensembl_gene_coordinate_dict[chrmsm][gene_putative_index + 1]]) is not list:
                        posterior_gene = self.gene_start_dict[chrmsm][self.ensembl_gene_coordinate_dict[chrmsm][gene_putative_index + 1]]
                        ens_gene_post = self.gene_id_exon_tuple_dict[posterior_gene]
                    else:
                        posterior_gene = None
                        ens_gene_post = None

                    actual_gene = get_actual_gene(self.ucsc_genes[genes],
                                                  ens_prior = ens_gene_prior, ens_post = ens_gene_post,
                                                  pre_gene = prior_gene, post_gene = posterior_gene)


            if actual_gene is not False:
                final_dict[actual_gene].update(self.ucsc_genes[genes])
            else:
                ucsc_transcripts_to_be_ignored.append(genes)




        return final_dict, ucsc_transcripts_to_be_ignored


    def make_new_boundaries(self):

        new_dict, ucsc_genes_to_be_ignored = self.combine_ucsc_specific_genes_to_ensembl_gene_dicts()
        gene_start_chr_dict_list = {}
        gene_start_gene_id_dict = {}
        temp_dict = {}

        for gene, transcript_dict in new_dict.items():
            try:
                chrmsm = self.ens_gene_id_chr_dict[gene]
            except KeyError:
                chrmsm = self.ucsc_gene_id_chr_dict[gene]
            start_list = []
            stop_list = []
            for transcript, exon_list in transcript_dict.items():
                for exons_tuple in exon_list:
                    start_list.append(int(exons_tuple[0]))
                    stop_list.append(int(exons_tuple[1]))
            left_boundary = min(start_list)
            right_boudary = max(stop_list)

            if not chrmsm in gene_start_chr_dict_list.keys():
                gene_start_chr_dict_list[chrmsm] = []

            if not chrmsm in gene_start_gene_id_dict.keys():
                gene_start_gene_id_dict[chrmsm] = {}

            try:
                if self.gene_id_direction_dict[gene] == "+":
                    gene_start_chr_dict_list[chrmsm].append(left_boundary)
                    gene_start_gene_id_dict[chrmsm][left_boundary] = gene
                else:
                    gene_start_chr_dict_list[chrmsm].append(right_boudary)
                    gene_start_gene_id_dict[chrmsm][right_boudary] = gene

            except KeyError:
                if self.ucsc_gene_direction_dict[gene] == "+":
                    gene_start_chr_dict_list[chrmsm].append(left_boundary)
                    gene_start_gene_id_dict[chrmsm][left_boundary] = gene
                else:
                    gene_start_chr_dict_list[chrmsm].append(right_boudary)
                    gene_start_gene_id_dict[chrmsm][right_boudary] = gene
        for chrm in gene_start_chr_dict_list.keys():
            relist_set_chr = list(set(gene_start_chr_dict_list[chrm]))
            temp_dict[chrm] = relist_set_chr

        gene_start_chr_dict_list = temp_dict

        for chr in gene_start_chr_dict_list.keys():
            gene_start_chr_dict_list[chr].sort()

        return gene_start_gene_id_dict, gene_start_chr_dict_list


    def tester(self):

        """
        One of the worst code. But in the interest of the correctness, I write the needful.
        Here I check if the miso coordinates are present. Sometimes, the coordinates are up or downstream of the
        putative gene and may fall into the downstream or the upstream gene.
        :return:
        """


        for chrmsm, gene_dict in self.miso_genedict.items():
            for miso_genes, isoforms in gene_dict.items():
                search_coordinate = self.gene_coordinate_search[miso_genes]
                gene_putative_location = gene_search(search_coordinate, self.gene_start_chr_dict_list[chrmsm],
                                                     self.miso_gene_direction_dict[miso_genes])
                #print(self.gene_start_gene_id_dict[chrmsm][gene_putative_location])
                gene_start_index = self.gene_start_chr_dict_list[chrmsm].index(gene_putative_location)
                try:
                    upstream_gene_coordinate = self.gene_start_chr_dict_list[chrmsm][gene_start_index-1]
                except IndexError:
                    upstream_gene_coordinate = gene_putative_location
                try:
                    downstream_gene_coordinate = self.gene_start_chr_dict_list[chrmsm][gene_start_index+1]
                except IndexError:
                    downstream_gene_coordinate = gene_putative_location
                gene, upstream_gene, downstream_gene = self.gene_start_gene_id_dict[chrmsm][gene_putative_location], \
                                                       self.gene_start_gene_id_dict[chrmsm][upstream_gene_coordinate], \
                                                       self.gene_start_gene_id_dict[chrmsm][downstream_gene_coordinate]
                gene_dict_pointer = self.final_dict[gene]
                upstream_gene_dict_pointer = self.final_dict[upstream_gene]
                downstream_gene_dict_pointer = self.final_dict[downstream_gene]
                for event in self.miso_geneid_dict.keys():
                    try:
                        event_dict = self.miso_geneid_dict[event][miso_genes]
                        self.identify_and_return_checked_genes(event, event_dict, gene, upstream_gene, downstream_gene,
                                                               gene_dict_pointer, upstream_gene_dict_pointer,
                                                               downstream_gene_dict_pointer)
                    except KeyError:
                        pass
        actual_gene_dict = self.gene_dict_with_MISO_features

        return actual_gene_dict


    def add_MISO_info(self, event, gene, transcript, coordinate = None):

        if not gene in self.gene_dict_with_MISO_features.keys():
            self.gene_dict_with_MISO_features[gene] = {}
            if not transcript in self.gene_dict_with_MISO_features[gene]:
                self.gene_dict_with_MISO_features[gene][transcript] = {}
                if not event in self.gene_dict_with_MISO_features[gene][transcript].keys():
                    self.gene_dict_with_MISO_features[gene][transcript][event] = []
                    if coordinate:
                        self.gene_dict_with_MISO_features[gene][transcript][event].append(coordinate)
                    else:
                        self.gene_dict_with_MISO_features[gene][transcript][event].append("Nothing")
                else:
                    if coordinate:
                        self.gene_dict_with_MISO_features[gene][transcript][event].append(coordinate)
                    else:
                        self.gene_dict_with_MISO_features[gene][transcript][event].append("Nothing")
            else:
                if not event in self.gene_dict_with_MISO_features[gene][transcript].keys():
                    self.gene_dict_with_MISO_features[gene][transcript][event] = []
                    if coordinate:
                        self.gene_dict_with_MISO_features[gene][transcript][event].append(coordinate)
                    else:
                        self.gene_dict_with_MISO_features[gene][transcript][event].append("Nothing")
                else:
                    if coordinate:
                        self.gene_dict_with_MISO_features[gene][transcript][event].append(coordinate)
                    else:
                        self.gene_dict_with_MISO_features[gene][transcript][event].append("Nothing")

        else:
            if not transcript in self.gene_dict_with_MISO_features[gene]:
                self.gene_dict_with_MISO_features[gene][transcript] = {}
                if not event in self.gene_dict_with_MISO_features[gene][transcript].keys():
                    self.gene_dict_with_MISO_features[gene][transcript][event] = []
                    if coordinate:
                        self.gene_dict_with_MISO_features[gene][transcript][event].append(coordinate)
                    else:
                        self.gene_dict_with_MISO_features[gene][transcript][event].append("Nothing")
                else:
                    if coordinate:
                        self.gene_dict_with_MISO_features[gene][transcript][event].append(coordinate)
                    else:
                        self.gene_dict_with_MISO_features[gene][transcript][event].append("Nothing")
            else:
                if not event in self.gene_dict_with_MISO_features[gene][transcript].keys():
                    self.gene_dict_with_MISO_features[gene][transcript][event] = []
                    if coordinate:
                        self.gene_dict_with_MISO_features[gene][transcript][event].append(coordinate)
                    else:
                        self.gene_dict_with_MISO_features[gene][transcript][event].append("Nothing")
                else:
                    if coordinate:
                        self.gene_dict_with_MISO_features[gene][transcript][event].append(coordinate)
                    else:
                        self.gene_dict_with_MISO_features[gene][transcript][event].append("Nothing")


    def identify_and_return_checked_genes(self, event, event_dict, gene, upstream_gene, downstream_gene,
                                          gene_dict_pointer, upstream_gene_dict_pointer, downstream_gene_dict_pointer):

        """
        USES RECURSION
        check if the important coordinates are present in the gene/transcripts;
        if not test if the upstream or the downstream coordinates are present in the gene identified.
        If none of the above works then test whether the coordinates are in the upstream or in the downstream genes.
        Even if nothing works for the above conditions, discard the MISO annotated coordinates.
        :param event:
        :param event_dict:
        :param gene_dict_pointer:
        :param upstream_gene_dict_pointer:
        :param downstream_gene_dict_pointer:
        :return:
        """
        i = 0
        not_found_dict = {}

        if event != "RI":
            for coordinates in event_dict["important_coordinates"]:

                if event == "A3SS" or event == "A5SS":
                    first, second = coordinates
                    for transcript, exon_list in gene_dict_pointer.items():
                        for exon in exon_list:
                            exon_start, exon_stop = exon
                            if int(first) == int(exon_start) or int(first) == int(exon_stop):
                                self.add_MISO_info(event, gene, transcript, coordinate = first)
                                i += 1
                            if int(second) == int(exon_start) or int(second) == int(exon_stop):
                                self.add_MISO_info(event, gene, transcript, coordinate = second)
                                i += 1

                    if i == 0:
                        if check_other_events_in_dict(event_dict, gene_dict_pointer):
                            self.add_MISO_info(event, gene, transcript)
                        elif check_other_events_in_dict(event_dict, upstream_gene_dict_pointer):
                            gene = upstream_gene
                            gene_dict_pointer = upstream_gene_dict_pointer
                            self.identify_and_return_checked_genes(event, event_dict, gene, upstream_gene,
                                                                   downstream_gene, gene_dict_pointer,
                                                                   upstream_gene_dict_pointer, downstream_gene_dict_pointer)
                        elif check_other_events_in_dict(event_dict, downstream_gene_dict_pointer):
                            gene = downstream_gene
                            gene_dict_pointer = downstream_gene_dict_pointer
                            self.identify_and_return_checked_genes(event, event_dict, gene, upstream_gene, downstream_gene,
                                                                   gene_dict_pointer, upstream_gene_dict_pointer,
                                                                   downstream_gene_dict_pointer)
                        else:
                            pass


                if event == "SE":
                    first, second = coordinates
                    for transcript, exon_list in gene_dict_pointer.items():
                        for exon in exon_list:
                            exon_start, exon_stop = exon
                            if int(first) == int(exon_start) and int(second) == int(exon_stop):
                                i += 1
                                self.add_MISO_info(event, gene, transcript, coordinate = exon)
                    if i == 0:
                        if check_other_events_in_dict(event_dict, gene_dict_pointer):
                            self.add_MISO_info(event, gene, transcript)
                        elif check_other_events_in_dict(event_dict, upstream_gene_dict_pointer):
                            gene = upstream_gene
                            gene_dict_pointer = upstream_gene_dict_pointer
                            self.identify_and_return_checked_genes(event, event_dict, gene, upstream_gene, downstream_gene,
                                                                   gene_dict_pointer, upstream_gene_dict_pointer,
                                                                   downstream_gene_dict_pointer)
                        elif check_other_events_in_dict(event_dict, downstream_gene_dict_pointer):
                            gene = downstream_gene
                            gene_dict_pointer = downstream_gene_dict_pointer
                            self.identify_and_return_checked_genes(event, event_dict, gene, upstream_gene, downstream_gene,
                                                                   gene_dict_pointer, upstream_gene_dict_pointer,
                                                                   downstream_gene_dict_pointer)
                        else:
                            pass

                if event == "MXE":
                    for coordinates in event_dict["important_coordinates"]:
                        for transcript, exon_list in gene_dict_pointer.items():
                            if coordinates in exon_list:
                                self.add_MISO_info(event, gene, transcript, coordinate = coordinates)
                                i += 1

                    if i == 0:
                        if check_other_events_in_dict(event_dict, gene_dict_pointer):
                            self.add_MISO_info(event, gene, transcript)
                        elif check_other_events_in_dict(event_dict, upstream_gene_dict_pointer):
                            gene = upstream_gene
                            gene_dict_pointer = upstream_gene_dict_pointer
                            self.identify_and_return_checked_genes(event, event_dict, gene, upstream_gene, downstream_gene,
                                                                   gene_dict_pointer, upstream_gene_dict_pointer,
                                                                   downstream_gene_dict_pointer)
                        elif check_other_events_in_dict(event_dict, downstream_gene_dict_pointer):
                            gene = downstream_gene
                            gene_dict_pointer = downstream_gene_dict_pointer
                            self.identify_and_return_checked_genes(event, event_dict, gene, upstream_gene, downstream_gene,
                                                                   gene_dict_pointer, upstream_gene_dict_pointer,
                                                                   downstream_gene_dict_pointer)
                        else:
                            pass


        if event == "RI":

            if not self.check_presence_of_RI_event(event, event_dict, gene, gene_dict_pointer):
                if not self.check_presence_of_RI_event(event, event_dict, upstream_gene, upstream_gene_dict_pointer):
                    if not self.check_presence_of_RI_event(event, event_dict,
                                                           downstream_gene, downstream_gene_dict_pointer):
                        pass


    def check_presence_of_RI_event(self, event, event_dict, gene, gene_dict_pointer):

        found = False
        ri_exon, non_ri_exon = check_presence_of_ri_intron(event_dict['important_coordinates'], event_dict['before'])
        additional_non_ri_exon = event_dict["after"][0]
        ri_exn1, ri_exn2 = ri_exon
        non_ri_exn1, non_ri_exn2 = non_ri_exon
        additional_non_ri_exn1, additional_non_ri_exn2 = additional_non_ri_exon

        for transcript, exon_list in gene_dict_pointer.items():
            non_ri_exon_list = []
            for exon in exon_list:
                exon1, exon2 = exon
                if int(ri_exn1) == int(exon1) and int(ri_exn2) == int(exon2):
                    self.add_MISO_info(event, gene, transcript, coordinate=exon)
                    found = True
                elif int(non_ri_exn1) == int(exon1) and int(non_ri_exn2) == int(exon2):
                    if not exon in non_ri_exon_list:
                        non_ri_exon_list.append(exon)
                        found = True
                else:
                    if int(additional_non_ri_exn1) == int(exon1) and int(additional_non_ri_exn2) == int(exon2):
                        if not exon in non_ri_exon_list:
                            non_ri_exon_list.append(exon)
                            found = True
        if found:
            return True
        else:
            return False

    def squeeze_unique_events(self):

        """
        Here we reduce the list of coordinates as they are repetitive in some case, especially in MXE events where the
        events work at the same coordinate but by diffferent mxe events.
        :return:
        """
        annotation_final_dict = {}
        for genes, transcript_dict in self.final_miso_dict.items():
            for transcript, event_dict in transcript_dict.items():
                for event, coordinate_list in event_dict.items():
                    if not genes in annotation_final_dict.keys():
                        annotation_final_dict[genes] = {}
                        annotation_final_dict[genes][transcript]={}
                        annotation_final_dict[genes][transcript][event] = list(set(coordinate_list))
                    else:
                        if not transcript in annotation_final_dict[genes].keys():
                            annotation_final_dict[genes][transcript] = {}
                            annotation_final_dict[genes][transcript][event] = list(set(coordinate_list))
                        else:
                            annotation_final_dict[genes][transcript][event] = list(set(coordinate_list))
        return annotation_final_dict


    def determine_start_stop_gene_boundaries_for_gtf(self):

        gene_boundaries_dict = {}
        transcript_boundaries_dict = {}

        for gene, transcript_dict in self.final_dict.items():
            for transcript, exon_list in transcript_dict.items():
                exon_start_list = []
                exon_stop_list = []
                for exon in exon_list:
                    exon_start_list.append(int(exon[0]))
                    exon_stop_list.append(int(exon[1]))

                start = min(exon_start_list)
                stop = max(exon_stop_list)
                transcript_boundaries_dict[transcript] = {"start":start, "stop":stop}
                if not gene in gene_boundaries_dict.keys():
                    gene_boundaries_dict[gene] = {}
                    gene_boundaries_dict[gene]["start"] = start
                    gene_boundaries_dict[gene]["stop"] = stop
                else:
                    if int(gene_boundaries_dict[gene]["start"]) > int(start):
                        gene_boundaries_dict[gene]["start"] = start
                    if int(gene_boundaries_dict[gene]["stop"]) < int(stop):
                        gene_boundaries_dict[gene]["stop"] = stop

        return gene_boundaries_dict, transcript_boundaries_dict


    def merge_descriptions_for_final_gtf_bak(self):

        fh = open(self.merged_fn, "w")

        chrmsm_list = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13",
                       "chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
        for chrmsm in chrmsm_list:
            for cd in self.gene_start_chr_dict_list[chrmsm]:
                try:
                    gene = self.gene_start_gene_id_dict[chrmsm][cd]
                except KeyError:
                    gene = None
                gene_direction = self.gene_direction[gene]
                gene_boundaries = self.gene_boundaries[gene]
                gene_start = gene_boundaries["start"]
                gene_stop = gene_boundaries["stop"]
                try:
                    gene_miso_dict_annotation = self.final_miso_annotation_dict[gene]
                except KeyError:
                    gene_miso_dict_annotation = None

                make_gene_annotation_for_writing(fh, chrmsm, gene, gene_direction, begin=gene_start,
                                                 end=gene_stop, annot=gene_miso_dict_annotation, feature="gene")
                for trscript, exn_list in self.final_dict[gene].items():
                    transcript_boundaries = self.transcript_boundaries[trscript]
                    transcript_start = transcript_boundaries["start"]
                    transcript_stop = transcript_boundaries["stop"]

                    if gene_miso_dict_annotation != None:
                        if trscript in gene_miso_dict_annotation.keys():
                            make_gene_annotation_for_writing(fh, chrmsm, gene, gene_direction, begin = transcript_start,
                                                             end = transcript_stop, transcript=trscript,
                                                             annot=gene_miso_dict_annotation, exon_list=exn_list,
                                                             feature="transcript")
                        else:
                            make_gene_annotation_for_writing(fh, chrmsm, gene, gene_direction, begin = transcript_start,
                                                             end = transcript_stop, transcript=trscript,
                                                             exon_list=exn_list, feature="transcript")
                    else:
                        make_gene_annotation_for_writing(fh, chrmsm, gene, gene_direction, begin = transcript_start,
                                                         end = transcript_stop, transcript=trscript,
                                                         exon_list=exn_list, feature="transcript")


        fh.close()




if __name__ == "__main__":
    data = CheckOverlaps(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    data.merge_descriptions_for_final_gtf_bak()
