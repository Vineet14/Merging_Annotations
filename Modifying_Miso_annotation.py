# With routine upgrades to the gene definitions there has been an out of syc nomenclature of the MISO annotations.
# Here I remove those MISO annotations that are absent from the merged GTF file. The output is going to be synced merged
# GTF as well as updated MISO files.

# author: Vineet Sharma
# Dated: 02/28/2020

from pprint import pprint as pp

def sorting_tuple(coordinate_tuple_list):

    new_list = []

    for indiv_tup in coordinate_tuple_list:
        a1, b1 = indiv_tup
        if int(a1) > int(b1):
            a1, b1 = b1, a1
        new_list.append((a1, b1))

    return new_list

def generating_tuples(string_coordinate_tuple):
    return_list = []
    tuple_strg = string_coordinate_tuple[1:-1]
    tuple_strg_list = tuple_strg.split("), (")
    if len(tuple_strg_list) > 1:
        first_element = tuple_strg_list[0][1:]
        first_element_list = first_element.split(", ")
        first_element_new_tuple = first_element_list[0][1:-1], first_element_list[1][1:-1]
        return_list.append(first_element_new_tuple)
        last_element = tuple_strg_list[len(tuple_strg_list)-1][:-1]
        last_element_list = last_element.split(", ")
        last_element_new_tuple = last_element_list[0][1:-1], last_element_list[1][1:-1]
        for elements in tuple_strg_list[1:-1]:
            new_str = elements[1:-1]
            list_new_str = new_str.split(", ")
            new_tuple = list_new_str[0][:-1], list_new_str[1][1:]
            return_list.append(new_tuple)
        return_list.append(last_element_new_tuple)
    else:
        element = tuple_strg_list[0]
        new_element_list = element.split(", ")
        single_tuple = new_element_list[0][2:-1], new_element_list[1][1:-2]
        return_list.append(single_tuple)

    return return_list

def test_presence_of_all_events(splice_event, coordinates_to_check, miso_geneid_dict, transcript_event_dict,
                                to_check=None, wherein=None, gn=None, cdn=None):

    counter_set = set()
    if splice_event in ["A3SS", "A5SS"]:
        for coordinate in coordinates_to_check:
            important_coordinate1, important_coordinate2 = coordinate
            for transcript, event_dict in transcript_event_dict.items():
                if splice_event in event_dict.keys():
                    if important_coordinate1 in event_dict[splice_event]:
                        counter_set.add(important_coordinate1)
                    if important_coordinate2 in event_dict[splice_event]:
                        counter_set.add(important_coordinate2)
        if len(counter_set) == 2:
            return True
        else:
            return False
    if splice_event == "SE":
        for coordinate in coordinates_to_check:
            for transcript, exon_list in wherein.items():
                if str(coordinate) not in exon_list:
                    if str(to_check["before"][0]) in exon_list:
                        if str(to_check["after"][0]) in exon_list:
                            return True
        return False
    if splice_event == "MXE":
        mxe_events = set()
        for coordinate in coordinates_to_check:
            for transcript, exon_list in wherein.items():
                if str(coordinate) in exon_list:
                    for upstream_coordinate in to_check["before"]:
                        if str(upstream_coordinate) in exon_list:
                            for downstream_coordinate in to_check["after"]:
                                if str(downstream_coordinate) in exon_list:
                                    mxe_events.add(str(coordinate))
        if len(list(mxe_events)) > 1:
            return True
        else:
            return False
    if splice_event == "RI":
        ri_events = set()
        for trns, exn_list in wherein.items():
            non_ri_cntr = 0
            for exn in exn_list:
                if cdn in miso_geneid_dict["A"]:
                    ri_events.add("A")
            for non_ri_exns in miso_geneid_dict["B"]:
                if str(non_ri_exns) in exn_list:
                    non_ri_cntr += 1
            if non_ri_cntr == 2:
                ri_events.add("B")

        if len(ri_events) > 1:
            return True
        else:
            return False


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

class ReadMergedFile(ReadGTFFile):

    def __init__(self, merged_fn):
        self.merged_fn = merged_fn
        ReadGTFFile.__init__(self, self.merged_fn)
        self.having_events_dict, self.gene_events_coordinate_dict, \
        self.gene_chromosome_dict, self.gene_transcript_event_dict, \
        self.gene_transcript_exon_dict = self.get_get_gene_nomeclature()

    def get_get_gene_nomeclature(self):
        no_events_dict = {}
        having_events_dict = {}
        gene_chromosome_dict = {}
        gene_events_coordinate_dict = {}
        gene_transcript_event_dict = {}
        gene_transcript_exon_dict = {}

        for lines in self.line_list:
            if lines[2] == "gene":
                chr = lines[0]
                description_data = lines[8].split(";")
                gene = description_data[0].split(" ")[1]
                gene_transcript_event_dict[gene] = {}
                gene_chromosome_dict[gene] = chr
                gene_events_coordinate_dict[gene] = {}
                gene_transcript_exon_dict[gene] = {}
                if description_data[1].startswith("events"):
                    events = description_data[1].split(" ")
                    gene_events = events[1].split("|")
                    gene_events[len(gene_events)-1] = gene_events[len(gene_events)-1][:-2]
                    having_events_dict.update({ gene: gene_events})
            if lines[2] == "transcript":
                tr_description_data = lines[8].split("; ")
                for transcript_data in tr_description_data:
                    if transcript_data.startswith("transcript_id"):
                        transcript = transcript_data.rstrip()[:-1].split(" ")[1]
                        gene_transcript_event_dict[gene][transcript] = {}
                        gene_transcript_exon_dict[gene][transcript] = []
                    if transcript_data.startswith("events"):
                        transcript_events = transcript_data.split("events_in_transcripts ")[1].rstrip()[:-1]
                        for tr_events in transcript_events.split("|"):
                            tr_event = tr_events.split(":")
                            name_event = tr_event[0]
                            coordinates = tr_event[1]
                            for coordinated in coordinates.split("#"):
                                coordinate = coordinated.lstrip().rstrip()
                                if name_event not in gene_transcript_event_dict[gene][transcript].keys():
                                    gene_transcript_event_dict[gene][transcript][name_event] = []
                                    gene_transcript_event_dict[gene][transcript][name_event].append(str(coordinate))
                                else:
                                    gene_transcript_event_dict[gene][transcript][name_event].append(str(coordinate))

                                if name_event not in gene_events_coordinate_dict[gene].keys():
                                    gene_events_coordinate_dict[gene][name_event] = set()
                                    gene_events_coordinate_dict[gene][name_event].add(coordinate.lstrip().rstrip())
                                else:
                                    gene_events_coordinate_dict[gene][name_event].add(coordinate.lstrip().strip())
            if lines[2] == "exon":
                exon_start = lines[3]
                exon_stop = lines[4]
                exon_tup = exon_start, exon_stop
                exon_str = str(exon_tup)
                gene_transcript_exon_dict[gene][transcript].append(exon_str)

        return having_events_dict, gene_events_coordinate_dict, gene_chromosome_dict, \
               gene_transcript_event_dict, gene_transcript_exon_dict


class SplicingMISOFeatures(ReadGTFFile):

    """
    Parses MISO annotations
    """

    def __init__(self, miso_fn):
        self.miso_fn = miso_fn
        self.miso_geneid_dict = {}
        ReadGTFFile.__init__(self, self.miso_fn)
        self.miso_genedict, self.miso_gene_direction_dict, self.miso_gene_id_chr_dict, \
        self.gene_coordinate_search, self.gene_id_lines = self.make_dict_MI()

    def make_dict_MI(self):

        gene_dict = {}
        gene_direction_dict = {}
        gene_id_chromosome_dict = {}
        gene_coordinate_look_dict = {}
        gene_id_lines = {}

        for lines in self.line_list:
            exon_tuple = lines[3], lines[4]
            chromosome = lines[0]
            miso_class = lines[1]
            for description in lines[8].split(";"):
                if description.lstrip().startswith("gene_id"):
                    gene_id = description.lstrip().split(" ")[1]
                    if gene_id not in gene_id_lines.keys():
                        gene_id_lines[gene_id] = []
                        gene_id_lines[gene_id].append(lines)
                    else:
                        gene_id_lines[gene_id].append(lines)
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


        return gene_dict, gene_direction_dict, gene_id_chromosome_dict, gene_coordinate_look_dict, gene_id_lines

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
                self.miso_geneid_dict["RI"][gene_id]["important_coordinates"].append((str(coordinates[1]), str(coordinates[2])))
                self.miso_geneid_dict["RI"][gene_id]["before"].append((str(coordinates[0]), str(coordinates[1])))
                self.miso_geneid_dict["RI"][gene_id]["before"].append((str(coordinates[0]), str(coordinates[3])))
                self.miso_geneid_dict["RI"][gene_id]["after"].append((str(coordinates[2]), str(coordinates[3])))


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

class MergingMergedMISO(ReadMergedFile, SplicingMISOFeatures):

    def __init__(self, merged_fn, miso_fn):
        self.merged_fn = merged_fn
        self.miso_fn = miso_fn
        ReadMergedFile.__init__(self, self.merged_fn)
        SplicingMISOFeatures.__init__(self, self.miso_fn)
        self.miso_id_to_keep = self.start_test()

    def start_test(self):
        events_to_be_kept = set()
        i = 0

        for gene, event_dict in self.gene_events_coordinate_dict.items():
            chr = self.gene_chromosome_dict[gene]
            if len(event_dict) != 0:
                for event, coordinates in event_dict.items():
                    for miso_id, miso_coordinate_dict in self.miso_geneid_dict[event].items():
                        if self.miso_gene_id_chr_dict[miso_id] == chr:
                            for coordinate in coordinates:
                                if event in ["A3SS","A5SS"]:
                                    try:
                                        for miso_coordinates in miso_coordinate_dict["important_coordinates"]:
                                            imp_coordinate1, imp_coordinate2 = miso_coordinates
                                            if int(coordinate) in [int(imp_coordinate1), int(imp_coordinate2)]:
                                                if test_presence_of_all_events(event,
                                                                               miso_coordinate_dict["important_coordinates"],
                                                                               self.miso_genedict[chr][miso_id],
                                                                               self.gene_transcript_event_dict[gene]):
                                                    events_to_be_kept.add(miso_id)
                                    except ValueError:
                                        pass

                                #if event in ["MXE", "SE"]:
                                if event == "SE":
                                    if coordinate != "Nothing":
                                        for imp_coordinates in miso_coordinate_dict["important_coordinates"]:
                                            cordn_list = generating_tuples(coordinate)
                                            for indiv_coordinate in cordn_list:
                                                if str(indiv_coordinate) == str(imp_coordinates):
                                                    if test_presence_of_all_events(event,
                                                                                   miso_coordinate_dict["important_coordinates"],
                                                                                   self.miso_genedict[chr][miso_id],
                                                                                   self.gene_transcript_event_dict[gene],
                                                                                   to_check=miso_coordinate_dict,
                                                                                   wherein=self.gene_transcript_exon_dict[gene],
                                                                                   gn=gene):
                                                        events_to_be_kept.add(miso_id)

                                if event == "MXE":
                                    if coordinate != "Nothing":
                                        for imp_coordinates in miso_coordinate_dict["important_coordinates"]:
                                            cordn_list = generating_tuples(coordinate)
                                            for indiv_coordinate in cordn_list:
                                                if str(indiv_coordinate) == str(imp_coordinates):
                                                    if test_presence_of_all_events(event,
                                                                                   miso_coordinate_dict["important_coordinates"],
                                                                                   self.miso_genedict[chr][miso_id],
                                                                                   self.gene_transcript_event_dict[gene],
                                                                                   to_check=miso_coordinate_dict,
                                                                                   wherein=self.gene_transcript_exon_dict[gene],
                                                                                   gn=gene):
                                                        events_to_be_kept.add(miso_id)

                                if event == "RI":
                                    for miso_cdn in miso_coordinate_dict["before"]:
                                        cordn_list = generating_tuples(coordinate)
                                        for indiv_coordinate in cordn_list:
                                            if str(indiv_coordinate) == str(miso_cdn):
                                                if test_presence_of_all_events(event,
                                                                               miso_coordinate_dict,
                                                                               self.miso_genedict[chr][miso_id],
                                                                               self.gene_transcript_event_dict[gene],
                                                                               to_check=miso_coordinate_dict,
                                                                               wherein=self.gene_transcript_exon_dict[gene],
                                                                               gn=gene, cdn = indiv_coordinate):

                                                    events_to_be_kept.add(miso_id)

        return events_to_be_kept

    def write_new_miso_file(self):

        fh = open("/Users/vs804/Desktop/HG19_new_all_miso_UCSC_Refseqmerged.gtf", "w")

        for gene, lines in self.gene_id_lines.items():
            if gene in self.miso_id_to_keep:
                for line in lines:
                    fh.write("\t".join(line))

        fh.close()





if __name__ == "__main__":

    data = MergingMergedMISO("/Users/vs804/Desktop/HG19_UCSC_ref_merged.gtf",
                             "/Users/vs804/PycharmProjects/HumanMergingGTFs/GTFs/hg19_MISO/all_Miso.gtf")
    data.write_new_miso_file()
