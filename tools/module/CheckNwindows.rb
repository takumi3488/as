#!/usr/bin/env ruby

#
#  check multi-window or not
#    written by HF
#
require 'KSUtil'
include KSUtil
class CheckNwindows
  attr_reader :nwin
  def initialize(filename)
    filename = parse_prefix(filename) + '.dat' unless filename.nil?
    if File.exist?(filename)
      @filename = filename
    else
      print "===== Error ======\n"
      print "file:#{filename} is not found.\n"
      exit
    end
    parse
  end

  def parse
    nwindow = false
    @nwin = nil
    linecount = nil
    open(@filename, 'r') do |fin|
      fin.each_line do |line|
        if nwindow
          linecount += 1
        end
        if line.chomp.strip =~ /^^NWIN/
          linecount = 0
          nwindow = true
        end
        if nwindow && linecount == 1
          @nwin = line.strip.chomp.split[0].to_i
        end
      end
    end
  end
end
