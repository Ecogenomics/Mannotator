#!/usr/bin/env ruby 

# == Synopsis 
#   takes a list of EC or KO numbers and produces annotated KEGG pathways
# == Examples 
#     markKeggPathway foo.txt
#
#   Other examples:
#     markKeggPathway -q bar.doc
#     markKeggPathway --verbose foo.html
# == Usage 
#   markKeggPathway [options] source_file
#   For help use: markKeggPathway.rb -h
# == Options
#   -h, --help          Displays help message
#   -v, --version       Display the version, then exit
#   -q, --quiet         Output as little as possible, overrides verbose
#   -V, --verbose       Verbose output
#   -k, --ko_numbers    Input is a list of KO numbers. Default is to use EC numbers
#   -a, --allpath       Evaluate every path that an enzyme is contained in.
#                       The default is to use only the first pathway returned by 
#                       KEGG, which is the most specific for that enzyme
#
# == Author
#   Connor Skennerton
# == Copyright
#    
#    Copyright (c) 2011 Connor Skennerton 
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.



require 'rubygems'
require 'optparse' 
require 'rdoc/usage'
require 'ostruct'
require 'date'
require 'bio'

class App
    VERSION = '0.0.1'
    
    attr_reader :options
    
    def initialize(arguments, stdin)
        @arguments = arguments
        @stdin = stdin
        @objs = []
        @queryPathways = []
        
        # Set defaults
        @options = OpenStruct.new
        @options.verbose = false
        @options.quiet = false
        @options.ko = false
        @options.allpath = false
        @serv = Bio::KEGG::API.new

        # TO DO - add additional defaults
    end
    
    # Parse options, check arguments, then process the command
    def run
        
        if parsed_options? && arguments_valid? 
            
            puts "Start at #{DateTime.now}\
            \
            " if @options.verbose
            
            output_options if @options.verbose # [Optional]
            
            process_arguments            

            read_ec_file
            
            
            if @options.ko
                @objs.each do |enzyme|
                    get_kegg_pathways_ko(enzyme)
                end
            else
                @objs.each do |enzyme|
                    get_kegg_pathways_ec(enzyme)
                end
            end
            
            uniq_pathways = get_unique_pathways(@queryPathways)
            
            uniq_pathways.each do |queryPath|
                mark_enzymes(queryPath)
            end
            #process_command
            
            puts "\
            Finished at #{DateTime.now}" if @options.verbose
            
            else
            output_usage
        end
        
    end
    
    protected
    
    def parsed_options?
        
        # Specify options
        opts = OptionParser.new 
        opts.on('-V', '--version')    { output_version ; exit 0 }
        opts.on('-h', '--help')       { output_help }
        opts.on('-v', '--verbose')    { @options.verbose = true }  
        opts.on('-q', '--quiet')      { @options.quiet = true }
        opts.on('-k', '--ko_number')  { @options.ko = true}
        opts.on('-a', '--allpath')    { @options.allpath = true}
        # TO DO - add additional options
        
        opts.parse!(@arguments) rescue return false
        
        process_options
        true      
    end
    
    # Performs post-parse processing on options
    def process_options
        @options.verbose = false if @options.quiet
    end
    
    def output_options
        puts "Options:\
        "
        
        @options.marshal_dump.each do |name, val|        
            puts "  #{name} = #{val}"
        end
    end
    
    # True if required arguments were provided
    def arguments_valid?
        # TO DO - implement your real logic here
        true if @arguments.length >= 1 
    end
    
    # Setup the arguments
    def process_arguments
        # TO DO - place in local vars, etc
    end
    
    def output_help
        output_version
        RDoc::usage() #exits app
    end
    
    def output_usage
        RDoc::usage('usage') # gets usage from comments above
    end
    
    def output_version
        puts "#{File.basename(__FILE__)} version #{VERSION}"
    end
    
    def read_ec_file
        File.open(@arguments[0], 'r') do |infile|
            while line = infile.gets
                @objs<<line.chomp
            end
        end
    end
    
    def get_kegg_pathways_ec(enzyme)
        pathways = @serv.get_pathways_by_enzymes(enzyme)
        if @options.verbose 
            puts "#{enzyme} is a member of #{pathways.length} pathway(s)"
        end
        if @options.allpath
            pathways.each do |path|
                @queryPathways<<path
            end
        else
            @queryPathways<<pathways.shift
        end
    end
    
    def get_kegg_pathways_ko(enzyme)
        pathways = @serv.get_pathways_by_kos(enzyme)
        if @options.verbose 
            puts "#{enzyme} is a member of #{pathways.length} pathway(s)"
        end
        if @options.allpath
            pathways.each do |path|
                @queryPathways<<path
            end
        else
            @queryPathways<<pathways.shift
        end
    end
        
        
    def get_unique_pathways(pathways)
        return pathways.uniq
    end
    
    def print_kegg_pathway(url2, queryPath)
        @serv.save_image(url2, "#{queryPath}.gif")
    end
    
    def mark_enzymes(queryPath)
        #queryPathways.each do |queryPath| 
        puts "marking enzymes in #{queryPath}" if options.verbose
        url2 = @serv.mark_pathway_by_objects(queryPath, @objs)
        print_kegg_pathway(url2, queryPath)
    end
    
    
    def process_command
        # TO DO - do whatever this app does
        
        #process_standard_input # [Optional]
    end
    
    def process_standard_input
        input = @stdin.read      
        # TO DO - process input
        
        # [Optional]
        # @stdin.each do |line| 
        #  # TO DO - process each line
        #end
    end
end


# TO DO - Add your Modules, Classes, etc


# Create and run the application
app = App.new(ARGV, STDIN)
app.run
