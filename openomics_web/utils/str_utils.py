COUNT = "_count"


def make_trie(words):
    root = dict()
    for word in words:
        current_dict = root
        for letter in word:
            current_dict = current_dict.setdefault(letter, {})

        if COUNT not in current_dict:
            current_dict[COUNT] = 1
        else:
            current_dict[COUNT] += 1
    return root


def longest_common_prefix(strs):
    trie = make_trie(strs)
    lcp_branches = [[letter, ] for letter in trie.keys()]

    def traverse_trie(dictionary, prefix):
        if len(dictionary.keys()) == 1 and COUNT not in dictionary.keys():
            letter = list(dictionary.keys())[0]
            traverse_trie(dictionary[letter], prefix + letter)
        else:
            return prefix

    for branch in lcp_branches:
        print("branch[0]", branch[0])
        branch.append(traverse_trie(trie, branch[0]))

    return lcp_branches
