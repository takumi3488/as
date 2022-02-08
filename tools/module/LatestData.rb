#!/usr/bin/env ruby

#
# search newest file
#
class LatestData
  def initialize(prefix, suffix, dir)
    @prefix = prefix
    @suffix = suffix
    @dir = dir
    @target = dir + '/' + prefix + '*' + suffix
  end

  def latestfile
    `ls -ltr #{@target}`.split("\n")[-1].split[-1]
  end

  def search_index
    File.basename(latestfile, @suffix).tr(@prefix, '').to_i
  end
end
