#!/usr/bin/env python3
import sys
import gettax3
import tax_analysis_parser3


def main():
    if __name__ != "__main__":
        return

    if len(sys.argv) < 3:
        print("need <stat hits file> <excluded ids file>", file=sys.stderr)
        return

    hits_file = sys.argv[1]
    excluded_ids = sys.argv[2]
    process_hits(hits_file, excluded_ids)


def process_hits(hits_file, id_file):
    resolved_hits = dict()
    excluded_ids = dict()
    scores = dict()
    lineage_cache = {}

    with gettax3.connect(None, "./gettax.sqlite", None, None) as conn:
        with open(hits_file) as f:
            for line in f:
                parts = line.strip('n').split('\t')
                spot = (parts[0])
                hits = set(int(p.split('x')[0]) for p in parts[1:])
                tax_id = tax_analysis_parser3.deduce_tax_id(hits, lineage_cache, conn)
                resolved_hits[spot] = tax_id

        with open(id_file) as f:
            for tid in f:
                lineage = tax_analysis_parser3.get_lineage(tid, lineage_cache, conn)
                excluded_ids[tid.strip('\n')] = lineage

        for spot in resolved_hits.keys():
            test_id = spot[6:].split('.')[0]
            hit_id = resolved_hits[spot]
            lin_len = len(excluded_ids[test_id])
            species_len = lin_len - 2
            genus_len = lin_len - 3
            is_virus = True if 10239 in excluded_ids[test_id] else False
            if hit_id == excluded_ids[test_id][species_len]:
                score = "species"
            elif hit_id == excluded_ids[test_id][genus_len]:
                score = "genus"
            elif hit_id != 1 and hit_id in excluded_ids[test_id]:
                score = "ancestor"
            else:
                lineage = tax_analysis_parser3.get_lineage(hit_id, lineage_cache, conn)
                sib_len = len(lineage)
                if lineage[(sib_len - 2)] == excluded_ids[test_id][genus_len]:
                    score = "species"
                elif hit_id != 1 and [value for value in lineage if value in excluded_ids[test_id] and value != 1]:
                    score = "ancestor"
                elif hit_id == 1:
                    score = "false negative"
                else:
                    score = "false positive"
                    print(test_id, spot, hit_id, file=sys.stderr)

            if is_virus:
                scores[spot] = ["Virus", score]
            else:
                scores[spot] = ["Bacteria", score]

        v_species = [[i for i in scores[x]] for x in scores.keys()
                     if 'species' in scores[x] and 'Virus' in scores[x]]
        v_genus = [[i for i in scores[x]] for x in scores.keys()
                   if 'genus' in scores[x] and 'Virus' in scores[x]]
        v_ancestor = [[i for i in scores[x]] for x in scores.keys()
                      if 'ancestor' in scores[x] and 'Virus' in scores[x]]
        v_fp = [[i for i in scores[x]] for x in scores.keys()
                if 'false positive' in scores[x] and 'Virus' in scores[x]]
        v_fn = [[i for i in scores[x]] for x in scores.keys()
                if 'false negative' in scores[x] and 'Virus' in scores[x]]
        b_species = [[i for i in scores[x]] for x in scores.keys()
                     if 'species' in scores[x] and 'Bacteria' in scores[x]]
        b_genus = [[i for i in scores[x]] for x in scores.keys()
                   if 'genus' in scores[x] and 'Bacteria' in scores[x]]
        b_ancestor = [[i for i in scores[x]] for x in scores.keys()
                      if 'ancestor' in scores[x] and 'Bacteria' in scores[x]]
        b_fp = [[i for i in scores[x]] for x in scores.keys()
                if 'false positive' in scores[x] and 'Bacteria' in scores[x]]
        b_fn = [[i for i in scores[x]] for x in scores.keys()
                if 'false negative' in scores[x] and 'Bacteria' in scores[x]]

        print("Virus:species: ", len(v_species))
        print("Virus:genus: ", len(v_genus))
        print("Virus:ancestor: ", len(v_ancestor))
        print("Virus:false positive: ", len(v_fp))
        print("Virus:false negative: ", len(v_fn))
        print("Bacteria:species: ", len(b_species))
        print("Bacteria:genus: ", len(b_genus))
        print("Bacteria:ancestor: ", len(b_ancestor))
        print("Bacteria:false positive: ", len(b_fp))
        print("Bacteria:false negative: ", len(b_fn))


main()
