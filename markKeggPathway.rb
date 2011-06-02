#!/usr/bin/env ruby

require 'rubygems'
require 'bio'

 objs = []
 File.open(ARGV[0], 'r') do |infile|
 	while line = infile.gets
 		objs<<line.chomp
 	end
 end
File.close
serv = Bio::KEGG::API.new
# mark pathway 

queryPathways = []

objs.each do |enzyme|
	pathways = serv.get_pathways_by_enzymes(enzyme)
	print "#{enzyme} is a member of "
	puts pathways.length
	#queryPathways#pathways.uniq!
	pathways.each do |path|
		#puts path
		queryPathways<<path
	end
end
queryPathways.uniq!
puts "querying KEGG for pathway images.  This will take alot of time!"
queryPathways.each do |queryPath| 
	puts "marking enzymes in #{queryPath}"
	url2 = serv.mark_pathway_by_objects(queryPath, objs)
	serv.save_image(url2, "#{queryPath}.gif")
end
# # save the result images
# serv.save_image(url1, "marked_pathway.gif")
# serv.save_image(url2, "Test_colored_pathway.gif")