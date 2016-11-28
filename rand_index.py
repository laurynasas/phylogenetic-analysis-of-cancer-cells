from numpy import *
from scipy.misc import comb


def get_number_of_same_objects_in_lists(a, b):
    counter = 0
    # print "a",a
    # print "b",b
    for el in a:
        if el in b:
            counter += 1
    return counter


def find_similarity_matrix(or_data, new_data):
    k = len(or_data.keys())
    # print or_data.keys()
    # print new_data.keys()
    sim_matrix = zeros((k, k))
    # print " \t",new_data
    # print " \t",or_data
    for i in xrange(k):
        for j in xrange(k):
            # print i,j
            sim_matrix[i, j] = get_number_of_same_objects_in_lists(or_data.get(i + 1), new_data.get(j))
    return sim_matrix


def find_common_labelling(or_data, new_data):
    sim_matrix = find_similarity_matrix(or_data, new_data)
    # print "----"
    # print sim_matrix
    # print "----"
    label_dict = {}
    or_labels = or_data.keys()
    for i, label in enumerate(or_labels):
        row= sim_matrix[i]
        points_to = argmax(row)
        while points_to in label_dict.values():
            row[points_to] = -1
            points_to = argmax(row)
        label_dict[label] = points_to
    return label_dict


def combine_old_new(labeling, original, new):
    combined = {}
    for key in original.keys():
        for old_list in original[key]:
            if combined.get(key):
                combined[key] += [old_list]
            else:
                combined[key] = [old_list]

        news = {}
        if key == 2:
            for el in combined[2]:
                if str(el) in news.keys():
                    news[str(el)] += 1
                else:
                    news[str(el)] = 1
            print "halfway combined counts", news
        # print labeling
        # print labeling[key]
        # print new
        for new_list in new[labeling[key]]:
            if combined.get(key):
                combined[key] += [new_list]
            else:
                combined[key] = [new_list]


    return combined


def calculate_tp_fp(labeling, original, new):
    k = len(labeling.keys())
    combined = combine_old_new(labeling, original, new)

    product = 0
    for key in combined.keys():
        product += comb(len(combined[key]), 2)

    return product


def find_index(el, original):
    for i, item in enumerate(original):
        if item == el:
            return i
    return None


def get_true_postive_score(original, combined, k):
    score = 0
    all_pairs = combined[k]
    for el in original[k]:
        all_pairs.pop(find_index(el, all_pairs))
        score += 1
    for el in original[k]:

        if el in all_pairs:
            all_pairs.pop(find_index(el, all_pairs))
            score += 1

    return score


def calculate_tp(labeling, original, new):
    combined = combine_old_new(labeling, original, new)
    product = 0
    for i in (combined.keys()):
        score = get_true_postive_score(original, combined, i)
        print " ",score
        if score > 1:
            product += comb(score, 2)

    return product


def get_scoring_of_lost_el(labeling, original, new, k_cluster):
    combined = combine_old_new(labeling, original, new)
    if k_cluster == 2:
        print "Original", original[2]
        print "new", new[labeling[2]]

        news = {}
        for el in new[labeling[2]]:
            if str(el) in news.keys():
                news[str(el)] += 1
            else:
                news[str(el)] = 1
        print "new:",news

        print "combined", combined[2]

    news ={}
    for el in combined[2]:
        if str(el) in news.keys():
            news[str(el)] += 1
        else:
            news[str(el)] = 1
    rabd ={}
    for el_master in new.values():
        for el in el_master:
            if str(el) in rabd.keys():
                rabd[str(el)] += 1
            else:
                rabd[str(el)] = 1

    if k_cluster == 2:
        print "combined counts: ", news
        # print rabd
    scores = {}
    for key in combined.keys():
        all_val = combined[key]

        for el in all_val:
            # print "   ",k_cluster
            if el in original[k_cluster]:
                if k_cluster == 2:
                    print "   element in:", el
                if scores.get(key):
                    scores[key] += 1
                else:
                    scores[key] = 1
                if k_cluster == 2:
                    print "so score is:", scores[key]
    return scores


def calculate_fn(labeling, original, new):
    overal_sum = 0
    for key in original.keys():
        scoring_lost_el = get_scoring_of_lost_el(labeling, original, new, key)
        print "->", scoring_lost_el, len(original[key]) * 2, len(new[labeling[key]])*2
        product = 0
        for cluster in scoring_lost_el.keys():
            if scoring_lost_el[cluster] > 1:
                product += comb(scoring_lost_el[cluster], 2)
        overal_sum += comb(len(original[key]) * 2, 2) - product
    return overal_sum


def calculate_tn(original, tp_fp, fn):
    total = 0
    for key in original.keys():
        total += len(original[key])
    return comb(total * 2, 2) - tp_fp - fn


def calculate_fp(tp_fp, tp):
    return tp_fp - tp


def rand_index(original, new):
    labeling = find_common_labelling(original, new)
    # print labeling
    tp = calculate_tp(labeling, original, new)
    print tp
    tp_fp = calculate_tp_fp(labeling, original, new)
    fn = calculate_fn(labeling, original, new)
    print fn
    tn = calculate_tn(original, tp_fp, fn)
    print tn
    fp = calculate_fp(tp_fp, tp)
    print fp
    return (tp + tn) / (tp + fp + fn + tn)


if __name__ == '__main__':
    original = {1: [[1, 0, 1, 1], [1, 1, 1, 1]],
                2: [[0, 0, 0, 1], [0, 0, 0, 0]],
                3: [[0, 0, 1, 1]]}
    #
    new = {0: [[0, 0, 0, 0], [0, 0, 0, 1]],
           1: [[1, 0, 1, 1], [1, 1, 1, 1]],
           2: [[0, 0, 1, 1]]}
    # labeling = {1: 1, 2: 2, 3: 3}

    print rand_index(original, new)
    # pass
