#!/usr/bin/env ruby

require 'rubygems'
require 'json'
require 'byebug'
require 'tempfile'

pdata = IO.read("protdata.json")

d = JSON.parse(pdata)

all_fasta_strings = {}
d.each do |pfamid,species|
  cur_fasta = String.new()
  #byebug
  species.each do |org,seq|
    begin
      #seq["nameline"] =~ /OS=([\w\s\.]+)[A-Z]{2}=/
      #nm = $1[0..-2].gsub(/\s/, "_")
      seq["nameline"] =~ /OS=(\w+)/
      nm = $1
    rescue
      print "uh oh... something went wrong"
      byebug
      print ""
    end
    cur_fasta << ">" + nm + "\n"
    cur_fasta << seq["seq"] + "\n"
  end
  all_fasta_strings[pfamid.to_sym] = cur_fasta
end

# Do alignments
all_fasta_strings.each do |k,v|
  tfile = Tempfile.new(["#{k.to_s}_", '.fasta'])

  tfile.write(v)
  tfile.rewind
  tfile.close

  `./muscle3.8.31_i86darwin64 -in #{tfile.path} -out aligned/#{k.to_s}.fasta`  

  tfile.unlink
end

byebug
puts ""
