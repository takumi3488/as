#!/usr/bin/env ruby
#
#

class ExecSCFseries
  def initialize(configure)
    @min_iter = configure[:min_iter]
    @max_iter = configure[:max_iter]
    @package = configure[:package]
    print_date
    go
  end

  def print_date
    print "Date: #{Time.now}\n"
    unless @min_iter == '0'
      print "#{@package} #{@min_iter} -> #{@max_iter}\n"
    end
  end

  def go
    t0 = Time.now
    @min_iter.to_i.upto(@max_iter.to_i) do |index|
      print "-> start #{@package} ... #{index}\n"
      t1 = Time.now
      `#{@package} #{index}`
      t2 = Time.new
      print "     end #{@package} ... #{index} ( #{t2 - t1}sec.)\n"
    end
    t9 = Time.now
    print "Total #{t9 - t0} (sec.)\n"
  end
end
