#!/usr/local/bin/python
#import libraries
import os
import argparse
from argparse import RawTextHelpFormatter
from collections import defaultdict
from lib.functional_group import Functional_group

#description to be displayed in the command line
desc = ("Project from Max Kellner.")
#initilze command line argument parser
parser = argparse.ArgumentParser(description = desc)
#add arguments to the ArgumentParser
#parser for organism1
parser.add_argument('-o1', '--organism1', dest = 'org1',
                    metavar = '"Homo sapiens"',
                    type = str,
                    action ='store',
                    help = 'Type the taxonomic name of an organism',
                    required = True)
#parser for organism2
parser.add_argument('-o2', '--organism2',
                    dest='org2',
                    metavar = '"Pan troglodytes"',
                    type = str,
                    action='store',
                    help = 'Type the taxonomic name of an organism',
                    required = True)
#parser for organisms to exlude, takes a list
parser.add_argument('-x', '--exclude',
                    dest = 'exclude',
                    metavar = '"Mus musculus"',
                    type = str,
                    nargs='+',
                    action = 'store',
                    help = 'Type in a list of organism that should be exluded \
                            in the orthologue search. Type "all" for orthologues \
                            specific between the 2 species.')
#parser for the input dir path
parser.add_argument('-i', '--input',
                    dest = 'input',
                    metavar = '"path"',
                    type = str,
                    action ='store',
                    help = 'Name of the input folder. Needs to contain the \
                            the files meNOG.members.tsv, meNOG.annotations.tsv,\
                            eggnog4.functional_categories.txt & \
                            eggnog4.species_list.txt.',
                    required = True)
#parser for the output file path
parser.add_argument('-o', '--output',
                    dest = 'output',
                    metavar = '"path"',
                    type = str,
                    action = 'store',
                    help = 'Files will be written to that path',
                    required = True)
#mutually exclusive group that requires true value, to avoid export conflicts
export = parser.add_mutually_exclusive_group(required=True)
#parser flags controlling should be written to the file
export.add_argument('-p', '--proteins',
                    dest = 'proteins',
                    action = 'store_true',
                    help = 'Write protein ids to the output path')
export.add_argument('-d', '--description',
                    dest = 'description',
                    action = 'store_true',
                    help = 'Write the number of proteins associated with a \
                            specific functional description to output path')
export.add_argument('-e', '--extensive',
                    dest = 'extensive',
                    action = 'store_true',
                    help = 'Write the orthologues group id, the associated \
                            proteins and the functional description \
                            to a file in output path')

#some functions...
def process_annotations_file(f):
    """Wrapper for read_file() that is designed to process the
    meNOG.annotations.tsv file and return a dictionary in the format
    d={Functional group id, description}.
    """
    group_col = 1
    desc_col = 5
    d=dict()
    lines = read_file(f) # read the file
    for l in lines:
        cols = l.strip().split('\t')
        func_grou_id = cols[group_col]
        desc = cols[desc_col]
        d[func_grou_id] = desc
    return d

def process_species_list_file(f):
    """Wrapper for read_file() that is designed to process the meNOG.members.tsv
    file and return a dictionary in the format d={Species name,taxon id}.
    """
    #some integers mapping columns to value
    name_col = 0
    taxid_col = 1
    d = dict()
    lines = read_file(f) # read the file
    for line in lines: #loop over the lines
        cols = line.strip().split('\t') #remove trailing whitespaces and split
        name = cols[name_col]
        tax_id = cols[taxid_col]
        d[name] = tax_id #assign values to dictionary
    return d

def process_members_file(f):
    """Function that passes the members file on on the parser, processes it and
    returns a dictionary in the format {Taxon ID, set(Functional_group)}.
    """
    lines = read_file(f) #read file
    d = defaultdict(set) #allows to directly add to the set via one command
    fg_list = []
    for l in lines:
        l = l.strip() #strip trailing \n
        # create a functional_group object that processes and stores one entry
        fg = Functional_group(l)
        d = add_functional_group_to_dict(fg, d) #add it to the dict
    return d

def process_functional_categories_file(f):
    """Function that passes the functional categories files on the parser,
    processes it and returns a dictionary in the format
    {key letter, function string}.
    """
    lines = read_file(f) #read the lines
    d = dict()
    header = "" #some var to hold the super group header
    #loop over the lines
    for l in lines:
        l = l.strip() #remove trailin whitespaces
        #when line is lenght 0 a new block begins, reset the header and jump to
        #the next line
        if l == 0:
            header = ""
            continue
        #when l is upper case assign it to the header
        header = "" if l.isupper() else header
        #add the more specific keyword to to the dictionary
        if l.startswith("["): #check if line contains key letter value paar.
            l = l.translate(None, '[]') #remove the brakets
            key, value = l.split(' ',1) #split at first whitespace
            value = header + " - " + value #add abstract header to specific value
            d[key] = value #add to dictionary
    return d

def add_functional_group_to_dict(fg, d):
    """Function that adds a Functional_group object as value to a dictionary,
    with all the tax_ids it contains als key and returns the dictionary.
    """
    #loop the tax_ids of one functional group
    for key in fg.tax_ids:
        d[key].add(fg) #add to the set in the dictionary
    return d

def read_file(f):
    """Read in a file and return all lines as list.
    """
    if os.path.isfile(f): #check if file exists
        with open(f, 'r') as fin:
            try:
                ("Parsing: " + f)
                lines = fin.readlines() # read file content
                fin.close()
                return lines
            except IOError: #if parsing failes throw error
                print "Could not read file: ", f
    else: #if file does not exist throw an error
        raise IOError("File does not exist: ", f)

def get_orthologues(fg_dict, org1, org2, organisms_exclude = None,):
    """Function that returns the number of homologue groups between species,
    their protein ids and the functional categories. Works with 2 organisms to
    include and n to exlude.
    """
    #add the functional groups known in org1 and 2 to a list
    fg_include = []
    fg_include.append(fg_dict[org1])
    fg_include.append(fg_dict[org2])
    #intersect to keep only the common ones
    fg_include_common = set.intersection(*fg_include)

    #if there are organisms to exclude, kick out the respective functional goups
    fg_exclude = []
    fg_exclude_common = set()
    #check if there to exclude
    if organisms_exclude != None and len(organisms_exclude) > 0:
        for o in organisms_exclude:
            fg_exclude.append(fg_dict[o])
        #perform a union to keep all groups that contain >1 of the organisms
        fg_exclude_common = set.union(*fg_exclude)

    fgs = fg_include_common - fg_exclude_common #substract the groups to exclude
    #get the proteins of the common functional groups for org1
    prots_per_organism = defaultdict(set)

    prots_per_category = defaultdict(set)

    for fg in fgs:
        prots = fg.get_proteins_for_species(org1)
        prots_per_category[fg.funct_cat].update(prots)
        prots_per_organism[org1].update(prots)

    #the number of functional groups is the number of orthologues
    number_orthologues = len(fgs)
    proteins_org1 = prots_per_organism[org1]

    #print some output for transparency..
    print "Number of ortholog groups: ", number_orthologues
    print "Number of proteins: ", len(proteins_org1)
    print ("Number of different functional categories (including combinations): "
            + str(len(prots_per_category)))
    #return everything as a tuple
    return fgs, number_orthologues, proteins_org1, prots_per_category

def write_protein_ids(org1, org2, exclude,
                      number_orthologues, prots_per_organism,
                      f):
    """Function that takes both organisms, the orthologues search results and a
    file path to which it writes the output.
    """
    excluded_organisms = ""
    if len(exclude) != 0 and exclude != None:
        excluded_organisms = " ".join(exclude)

    header = ("# " + str(number_orthologues) +
              " orthologues of " + org1 + " in " + org2)
    header = header + " - excluded_organisms: " + excluded_organisms + "\n"
    header = header + "protein_ids"
    lines = []
    lines.append(header)

    sorted_prots = sorted(list(prots_per_organism))
    ids_to_write = '\n'.join(sorted_prots).rstrip()
    lines.append(ids_to_write)
    write_to_file(f,lines)

def write_functional_groups_counts(catgories_prots, categories_desc,f):
    """Function that takes dictionary mapping proteins to functional catagories
    id, a dictionary mapping functional catagories id to a string and a file
    output name.
    """
    d = dict()
    for key in catgories_prots:
        members = catgories_prots[key]
        numb_members = len(members)
        description = []
        #key can be combination of chars, each representing one category
        for c in key:
            description.append(categories_desc[c])
        description = ';'.join(description).strip(';')
        line = str(numb_members) + "\t" + description
        d[numb_members] = line

    lines = []
    for key in sorted(d):
        lines.append(d[key])
    write_to_file(f,lines)

def write_exensive(fgs, org1, org2, anno_dict,f):
    lines = []
    for fg in fgs:
        group_id = fg.group_id #get the orthologue group id
        #get the proteins for one organism, convert to list and sort
        prots_org1 = sorted(list(fg.get_proteins_for_species(org1)))
        #generate string where the proteins are delimited by ','
        #and remove the trailng ','
        prots_org1 = ','.join(prots_org1).strip(',')
        prots_org2 = sorted(list(fg.get_proteins_for_species(org2))) #see above
        prots_org2 = ','.join(prots_org2).strip(',') #see above
        annotation = anno_dict[group_id] #get the corresponding annotation

        #join everything to a line that should be written to file
        line = fg.group_id + '\t' + prots_org1 + '\t' + prots_org2 + '\t' + \
               annotation
        lines.append(line) #add to all lines
    write_to_file(f,lines) #write to file

def write_to_file(f,lines):
    """Function that takes a file path and a list as input. Will try to
    write the list as lines to the specified filed path.
    """
    with open(f, 'w') as fout:
        try:
            for l in lines:
                print l
                fout.write(l + '\n')
            fout.close()
        except IOError:
            print "Could not write to file: ", f

#main method
def main():
    """The main method to be excuted when the script is run. The main method
    will first parse command line arguments, then map organsism names to
    tax ids, then use these to filter the orthologues groups and then run the
    approriate commands for writing the expected output to a file.
    """
    #parse all the arguments given
    args = parser.parse_args()
    #assign the arguments to variables
    org1 = args.org1
    org2 = args.org2
    org_exlude = args.exclude

    f_in = args.input
    f_out = args.output

    store_prots = args.proteins
    store_desc = args.description
    store_extensive = args.extensive

    #define some file paths
    members_file = os.path.join(f_in, 'meNOG.members.tsv')
    annotations_file = os.path.join(f_in, 'meNOG.annotations.tsv')
    categories_file = os.path.join(f_in, 'eggnog4.functional_categories.txt')
    species_file = os.path.join(f_in, 'eggnog4.species_list.txt')

    #but do some checks before trying to parsing everything
    if org1 == org2:
        raise ValueError("Organism 1 and 2 can't be the same.")

    #parse the files and store the data in dictionaries
    organism_dict = process_species_list_file(species_file)
    fg_dict = process_members_file(members_file)

    #get the taxon id of the organisms
    tax_id_1 = organism_dict[org1]
    tax_id_2 = organism_dict[org2]

    #create a blacklist of organisms if supplied...
    organism_blacklist = set()
    if org_exlude is not None:
        if 'all' in org_exlude:
            organism_blacklist = set(organism_dict.values())
            organism_blacklist = organism_blacklist - set([tax_id_1, tax_id_2])
        else:
            for o in org_exlude:
                 if o not in organism_dict:
                     raise LookupError(o, "Not included in the database... ")
                 else:
                     organism_blacklist.add(organism_dict[o])

    #get the number of orthologues groups and the proteins for org1
    fgs, number_orthologues, proteins, category_dict = get_orthologues(fg_dict,
                                                                     tax_id_1,
                                                                     tax_id_2,
                                                                     organism_blacklist)
    #run the commands to get the respective meta data and write results to files
    if store_prots:
        write_protein_ids(org1, org2, org_exlude, number_orthologues, proteins,
                          f_out)
    if store_desc:
        category_mapping_dict = process_functional_categories_file(categories_file)
        write_functional_groups_counts(category_dict,category_mapping_dict,f_out)

    if store_extensive:
        anno_dict = process_annotations_file(annotations_file)
        write_exensive(fgs, tax_id_1, tax_id_2 , anno_dict, f_out)
#make code only excutable when run as program
if __name__ == "__main__": main()
