#!/usr/bin/env ruby

require 'optparse'
require 'optparse/time'
require 'ostruct'
require 'pp'
require 'set'

class OptionParse
	#
	# Return a structure describing the options.
	#
	def self.parse(args)
		# The options specified on the command line will be collected in *options*.
		# We set default values here.
		options = OpenStruct.new
		options.read_count = 1.0
		
		opts = OptionParser.new do |opts|
			opts.banner = "Usage: example.rb [options]"
		
			opts.separator ""
			opts.separator "Specific options:"
		
			# Mandatory argument.
			opts.on("-g", "--gff_file FILE",
				"") do |gff_file|
				options.gff_file_name = gff_file
			end
			opts.on("-b", "--bed_file FILE",
				"") do |bed_file|
				options.bed_file_name = bed_file
			end
		end
		
		opts.parse!(args)
		options
	end  # parse()

end  # class OptparseExample


class GffFeature
	@id
	@start_coord
	@end_coord
	@strand

	def initialize(id, start_coord, end_coord, strand)
		@id=id
		@start_coord=start_coord
		@end_coord=end_coord
		@strand=strand
	end

	def <=>(other)
		self.oriented_start <=> other.oriented_start
	end


	def oriented_start
		if @strand=="+"
			return @start_coord
		elsif @strand=="-"
			return @end_coord
		else
			return @start_coord
		end
	end

	#Accessors
	def id
		return @id
	end

	def start_coord
		return @start_coord
	end

	def end_coord
		return @end_coord
	end

	def strand
		return @strand
	end

end

def rangeOrderedSet(start_coord,end_coord,target_set)
	
	set_size = target_set.length
	lower = 0
	upper = set_size-1
	while lower + 1 != upper
		mid = ((lower+upper)/2).to_i
		if target_set[mid].oriented_start < start_coord
			lower=mid
		else
			upper=mid
		end
	end
	if target_set[lower].oriented_start > start_coord
		begin_range = lower
	else
		begin_range = upper
	end

	lower = 0
	upper = set_size-1
	while lower + 1 != upper
		mid = ((lower+upper)/2).to_i
		if target_set[mid].oriented_start < end_coord
			lower=mid
		else
			upper=mid
		end
	end
	if target_set[upper].oriented_start < end_coord
		end_range = upper
	else
		end_range = lower
	end


	if begin_range > end_range
		return Array.new
	end

	return target_set[begin_range..end_range]
end

seqids = Hash.new

options = OptionParse.parse(ARGV)

# Process GFF file
gff_file = File.open(options.gff_file_name)
gff_file.each_with_index do |line,i|
	if line.strip[0]!="#"
		segline = line.split("\t")
		seqid = segline[0].strip
		if seqids[seqid].nil?
			seqids[seqid] = Array.new
		end
		if segline[2]=="gene"
			feature = GffFeature.new(segline[8].split(";")[0].split("=")[1],segline[3].to_i,segline[4].to_i,segline[6])
			seqids[seqid] << feature
		end
	end
end
seqids.each_key {|key| seqids[key].sort! { |a,b| a <=> b }}

# Process BED file
bed_file = File.open(options.bed_file_name)
bed_file.each_with_index do |line,i|
	segline = line.split("\t")
	seqid = segline[0].strip
	start_coord = segline[1].strip.to_i
	end_coord = segline[2].strip.to_i
	mid_coord = ((start_coord+end_coord)/2).to_i
	neighbor_features = rangeOrderedSet(mid_coord-1000,mid_coord+1000,seqids[seqid])
	puts line

	location="INTERGENIC"
	tss_distance=1000
	nearest_feature=nil
	neighbor_features.each_with_index do |feature,i|
		if (feature.strand=="+" && (mid_coord - feature.start_coord) > -tss_distance && (mid_coord - feature.end_coord) <= 0) || (feature.strand=="-" && (mid_coord - feature.end_coord) < tss_distance && (mid_coord - feature.start_coord) >= 0)
			tss_distance = mid_coord - feature.oriented_start
			if feature.strand=="-"
				tss_distance = -1 * tss_distance
			end
			nearest_feature = feature
		end
		if mid_coord > feature.start_coord && mid_coord < feature.end_coord
			location="GENE_BODY"
		end
	end
	if neighbor_features.size == 0 || nearest_feature.nil?
		puts "NEIGHBOR_TSS_COUNT: " + neighbor_features.size.to_s
	else
		puts "NEIGHBOR_TSS_COUNT: " + neighbor_features.size.to_s + "\t" + location + "\t" + nearest_feature.strand  + "\t" + tss_distance.to_s
	end

	neighbor_features.each_with_index do |feature,i|
		puts feature.id
	end
end

#chr1_subset = rangeOrderedSet(6000,7000,seqids["chr1"])
#chr1_subset.each_with_index do |feature,i|
#	puts feature.id
#end

