require 'rubygems'
require 'tempfile'
require 'find'
require 'byebug'

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
