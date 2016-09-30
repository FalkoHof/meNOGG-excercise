from collections import defaultdict

class Functional_group:

  def __init__(self, line):
      #some integers mapping columns to value
      self.group_id_col = 1
      self.prot_count_col = 2
      self.species_count_col = 3
      self.funct_cat_col = 4
      self.tax_prot_ids_col = 5
      #split the line at tab
      self.cols = line.strip().split('\t')
      #assign some values based on the column index
      self.group_id = self.cols[self.group_id_col]
      self.prot_count = self.cols[self.prot_count_col]
      self.species_count = self.cols[self.species_count_col]
      self.funct_cat = self.cols[self.funct_cat_col]
      #split taxid.protein_ids string into single taxid.protein_ids
      self.tax_prot_ids = set(self.cols[self.tax_prot_ids_col].split(','))
      #split taxid.protein_ids into taxid,protein ids and get a dictionary
      #in the form of {tax_ids,set(protein_ids)}
      self.tax_ids, self.prot_ids, self.tax_prot_dict = \
      self.process_tax_prot_ids(self.tax_prot_ids)

  def process_tax_prot_ids(self, tax_prot_ids):
      """Function that takes a list of ["taxid1.protein_id1", ..] and splits
      the string it the first dot. Returns a tuple of set(tax_ids),
      set(protein_ids), dict{tax_id,set(protein_ids)}
      """
      #initilze some objects to hold data
      tax_prot_dict = defaultdict(set)
      tax_ids = set()
      prot_ids = set()
      #iterate over all functional group members and process them
      for i in tax_prot_ids:
          #split only at the first . as some ids contain dots
          tax_id, prot_id = i.split('.', 1)
          tax_ids.add(tax_id)
          prot_ids.add(prot_id)
          tax_prot_dict[tax_id].add(prot_id)
      return tax_ids, prot_ids, tax_prot_dict

  def get_proteins_for_species(self, tax_id):
      """Function that takes a taxon id and reuturns the proteins assosiated
      with that species as set. If the species is not present throw an error.
      """
      if not tax_id in self.tax_prot_dict:
          raise LookupError("Taxon id not present in " + self.group_id)
      else:
          return self.tax_prot_dict[tax_id]
