#!/usr/bin/env ruby

require 'rubygems'
require 'find'
require 'byebug'
require 'tempfile'
require 'csv'

#########################################################################
# Convert fasta files to ruby objects (AND SORT GENUSES ALPHABETICALLY) #
#########################################################################
# find fasta files
fasta_file_paths = []
Find.find('aligned') { |path| fasta_file_paths << path if path =~ /\/[A-z0-9]+\.fasta$/ }
reference_fasta_path = 'aligned/aligned_cyt_b.fasta'
fasta_file_paths.reject! { |i| i == reference_fasta_path }


# parse and sort reference peptides
ref_fasta_parts = IO.read(reference_fasta_path).split('>')[1..-1]
#ref_trimmed = ref_fasta_parts.map { |x| x.split(') [')[1].chomp() } #<== Relic from old ver
ref_trimmed = ref_fasta_parts
ref_seqs = ref_trimmed.map { |s| {:genus => s.lines.first.chomp(), :seq => s.split("\n")[1..-1].join.gsub("\n",'') } }
ref_seqs.sort_by! { |h| h[:genus] }

# parse and sort other families
def fasta_parse(fasta_string)
  fasta_parts = fasta_string.split('>')[1..-1]
  #no_pfamid = fasta_parts.map { |x| x.split(',')[1].chomp() } #<== Relic from old ver
  no_pfamid = fasta_parts

  seqs = no_pfamid.map { |s| {:genus => s.lines.first.chomp(), :seq => s.split("\n")[1..-1].join.gsub("\n",'') } }
  seqs.sort_by! { |h| h[:genus] } #makes sure to sort alphabetically, allowing us to build partitioned character matrices
  return seqs
end
prot_families = fasta_file_paths.map { |f| fasta_parse(IO.read(f)) }

# Get lists of samples size 11
combs = []
1000.times { |_| combs.push((1..49).to_a.sample(12).sort()) }


rf_results ||= "" # only initialize variable if not yet defined
rf_values ||= []

sample_num = 0
combs.each do |c|
  sample_num += 1

  #################################################
  # # Build tnt/hennig format files for tnt input #
  #################################################
  def make_tnt_input(seq_arr, c)
    seq_arr_subset = []
    c.each { |x| seq_arr_subset.push(seq_arr[x-1]) } # `x-1` since seq_arr initializes from 0
    #byebug
    seq_arr_length = seq_arr_subset[0][:seq].length()
    retstr = "xread\n#{seq_arr_length} 12\n"  #<== NOTE!!!!!! Must be set to number of genera!!!!!!
    begin
      seq_arr_subset.each { |sa| retstr << "#{sa[:genus]} #{sa[:seq]}\n" }
    rescue
      byebug
      puts ""
    end
    retstr << ";\n"
  end
  tnt_input_strings = prot_families.map { |pf| make_tnt_input(pf, c) }
  reftree_input_string = make_tnt_input(ref_seqs, c)

  ######################
  # # Run TNT analysis #
  ######################
  # Iterate over all protein families plus reference tree; make trees using aquickie TNT script, save to appropriate file name
  $i = 0
  until $i > 4 do
    name = fasta_file_paths[$i].split("/")[-1].split(".")[0]
    file = Tempfile.new(['tntinput','.tnt']) # creates tempfile with .tnt extension
    File.open(file.path, 'w') {|f| f.write(tnt_input_strings[$i]) }
    #byebug
    begin
      `./tnt.command mxram 4000, nstate prot,proc #{file.path},jdrquickie,zzz,`
      # do something with the tree that was output
      `mv jdrquickie.svg trees/#{name}.svg`
      `mv aquickie.tre trees/tree_files/#{name}.tre`
    ensure
      file.close
      file.unlink
    end
    $i += 1
  end
  # Make tree for reference
  name = "reference_cyt_b"
  #file = Tempfile.new(['tntinput','.tnt'])
  file = 'referencedata.tnt'
  File.open(file,'w') {|f| f.write(reftree_input_string) }
  begin
    `./tnt.command nstate prot, proc #{file},jdrquickie,zzz,`
    `mv jdrquickie.svg trees/ref_cyt_b.svg`
    `mv aquickie.tre trees/tree_files/ref_cyt_b.tre`
  ensure
    # file.close
    # file.unlink
  end

  # Get tree file names
  tree_file_paths = []
  Find.find('trees/tree_files') { |path| tree_file_paths << path if path =~ /\.tre$/ }

  rf_combs = tree_file_paths.combination(2).to_a
  rf_combs.reject! { |i| ! i.include?("trees/tree_files/ref_cyt_b.tre") }

  rf_combs.each do |comb|
    ntrees1 = (`wc -l #{comb[0]}`.to_i - 2)
    ntrees2 = (`wc -l #{comb[1]}`.to_i - 2)
    
    `./tnt.command nstates prot, p referencedata.tnt, p #{comb[0]}, p #{comb[1]}, rfdistances 0 1,zzz,`
    rf_results << "#{sample_num} #{comb[0].split('/')[-1].split('.')[0]} #{comb[1].split('/')[-1].split('.')[0]} #{IO.read('rflog.log')}"
    rf_values.push(IO.read('rflog.log').to_s)
  end

  # distance of top two trees in each set - UNNECESSARY
  # tree_file_paths.each do |tfp|
  #   `./tnt.command nstates prot, p referencedata.tnt, p #{tfp}, rfdistances 0 1,zzz`
  #   rf_results << "#{tfp.split('/')[-1].split('.')[0]} #{IO.read('rflog.log')}"
  # end

end

CSV.open("all_rf_vals.csv", 'wb') do |csv|
  rf_results.each_line { |l| csv << l.split(" ") }
end

puts 'All done!'
