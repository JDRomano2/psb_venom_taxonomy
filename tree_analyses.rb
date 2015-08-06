require 'rubygems'
require 'bio'
require 'json'
require 'find' # finds files
require 'byebug'
require 'tempfile'

#########################################################################
# Convert fasta files to ruby objects (AND SORT GENUSES ALPHABETICALLY) #
#########################################################################
# find fasta files
fasta_file_paths = []
Find.find('aligned') { |path| fasta_file_paths << path if path =~ /\/[A-z0-9]+.muscle\.fasta$/ }
reference_fasta_path = 'aligned/aligned_cyt_b.fasta'

byebug

# parse and sort reference peptides
ref_fasta_parts = IO.read(reference_fasta_path).split('>')[1..-1]
#ref_trimmed = ref_fasta_parts.map { |x| x.split(') [')[1].chomp() }
ref_trimmed = ref_fasta_parts
ref_seqs = ref_trimmed.map { |s| {:genus => s.lines.first.chomp().gsub("]",'').split(" ")[0], :seq => s.split("\n")[1..-1].join.gsub("\n",'') } }
ref_seqs.sort_by! { |h| h[:genus] }

# get random samples of 11 1000 times from 49 elements (# species in cyt_b matrices)
combs = []
10.times { |_| combs.push((1..49).to_a.sample(11)) }

# For 


byebug
puts ""


# parse and sort other families
def fasta_parse(fasta_string)
  fasta_parts = fasta_string.split('>')[1..-1]
  no_pfamid = fasta_parts.map { |x| x.split(',')[1].chomp() }
  seqs = no_pfamid.map { |s| {:genus => s.lines.first.chomp(), :seq => s.split("\n")[1..-1].join.gsub("\n",'') } }
  seqs.sort_by! { |h| h[:genus] } #makes sure to sort alphabetically, allowing us to build partitioned character matrices
  return seqs
end
prot_families = fasta_file_paths.map { |f| fasta_parse(IO.read(f)) }

#################################################
# # Build tnt/hennig format files for tnt input #
#################################################
def make_tnt_input(seq_arr)
  seq_arr_length = seq_arr[0][:seq].length()
  retstr = "xread\n#{seq_arr_length} 49\n"
  seq_arr.each { |sa| retstr << "#{sa[:genus]} #{sa[:seq]}\n" }
  retstr << ";\n"
end
tnt_input_strings = prot_families.map { |pf| make_tnt_input(pf) }
reftree_input_string = make_tnt_input(ref_seqs)

######################
# # Run TNT analysis #
######################
# Iterate over all protein families plus reference tree; make trees using aquickie TNT script, save to appropriate file name
$i = 0
# until $i > 3 do
#   name = fasta_file_paths[$i].split("/")[-1].split(".")[0]
#   file = Tempfile.new(['tntinput','.tnt']) # creates tempfile with .tnt extension
#   File.open(file.path, 'w') {|f| f.write(tnt_input_strings[$i]) }
#   begin
#     `sequences/aligned/tnt.command nstate prot,proc #{file.path},jdrquickie,zzz,`
#   # do something with the tree that was output
#     `mv jdrquickie.svg trees/#{name}.svg`
#     `mv tmp.tre trees/tree_files/#{name}.tre`
#   ensure
#     file.close
#     file.unlink
#   end
#   $i += 1
# end
# Make tree for reference
name = "reference_cyt_b"
#file = Tempfile.new(['tntinput','.tnt'])
file = 'referencedata.tnt'
File.open(file,'w') {|f| f.write(reftree_input_string) }
begin
  `./tnt.command nstate prot, proc #{file},jdrquickie,zzz,`
  `mv jdrquickie.svg trees/ref_cyt_b.svg`
  `mv tmp.tre trees/tree_files/ref_cyt_b.tre`
ensure
  # file.close
  # file.unlink
end


#####################################################################
# Compute incongruence length difference between all pairs of trees #
#####################################################################

#ref_block = fasta_parse(IO.read('sequences/aligned/aligned_cyt_b.fasta'))

def do_ild(block1, block2)

  total_length = block1[0][:seq].length() + block2[0][:seq].length()
  ildstring = ""
  ildstring << "nstates prot;\n"
  ildstring << "xread\n"
  ildstring << "#{total_length} 49\n"
  ildstring << "\n"
  ildstring << "&[prot]\n"
  block1.each { |b1| ildstring << "#{b1[:genus]} #{b1[:seq]}\n" }
  ildstring << "\n"
  ildstring << "&[prot]\n"
  block2.each { |b2| ildstring << "#{b2[:genus]} #{b2[:seq]}\n" }

  ildfile = Tempfile.new(['ildinput','.tnt'])
  File.open(ildfile.path,'w') { |f| f.write(ildstring) }
  `sequences/aligned/tnt.command mxram 2000, sect:slack 35, p #{ildfile.path}, ild, zzz`
  puts "p-value?"
  pval = gets.chomp()
  return pval

end
# combs = (0..3).to_a.combination(2).to_a
# res = []
# combs.each { |c| res << do_ild(prot_families[c[0]],prot_families[c[1]]) }
# res = []
# (0..3).to_a.each { |c| res << do_ild(ref_seqs, prot_families[c]) }




# Get tree file names
tree_file_paths = []
Find.find('trees/tree_files') { |path| tree_file_paths << path if path =~ /\.tre$/ }

rf_combs = tree_file_paths.combination(2).to_a
rf_results = ""
rf_combs.each do |comb|
  `./tnt.command nstates prot, p referencedata.tnt, p #{comb[0]}, p #{comb[1]}, rfdistances 0 2,zzz,`
  rf_results << "#{comb[0].split('/')[-1].split('.')[0]} #{comb[1].split('/')[-1].split('.')[0]} #{IO.read('rflog.log')}"
end

tree_file_paths.each do |tfp|
  `./tnt.command nstates prot, p referencedata.tnt, p #{tfp}, rfdistances 0 1,zzz`
  rf_results << "#{tfp.split('/')[-1].split('.')[0]} #{IO.read('rflog.log')}"
end




byebug
puts ''
