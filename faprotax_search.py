import re

def parse_faprotax(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    data = {}
    current_group = None

    for line in lines:
        # Skip comments and empty lines
        if line.startswith('#') or line.strip() == '':
            continue

        # Check if the line is a group header
        if not line.startswith('*'):
            current_group = line.strip()
            data[current_group] = []
        else:
            if current_group:
                data[current_group].append(line.strip())

    return data

def search_faprotax(data, search_term):
    result = {}
    for group, taxa in data.items():
        matched_taxa = []
        for taxon in taxa:
            taxon_part = taxon.split('#')[0].strip()
            if search_term.lower() in taxon_part.lower():
                matched_taxa.append(taxon)
        if matched_taxa:
            result[group] = matched_taxa

    return result

def print_search_results(results):
    for group, taxa in results.items():
        print(f"\n\n{group}")
        for taxon in taxa:
            print([i.strip() for i in taxon.split('# ')])

# file_path = 'FAPROTAX_1.2.10/FAPROTAX.txt'  # Replace with your file path
# data = parse_faprotax(file_path)
#
# search_term = 'Janthinobacterium'
# results = search_faprotax(data, search_term)
#
# # Print the results
# print_search_results(results)


# %%

def faprotax_search(search_term, file_path='FAPROTAX_1.2.10/FAPROTAX.txt'):
    data = parse_faprotax(file_path)
    return [i.split('\t', 1)[0] for i in list(search_faprotax(data, search_term).keys())]
