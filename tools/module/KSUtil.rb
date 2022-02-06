#!/usr/bin/env ruby
#
#
module KSUtil
  def parse_prefix(filename)
    prefix = filename
    ext_list = %w(dst.dat eld.dat esd.dat chg.dat fchg.dat pot.dat fpot.dat )
    ext_list.each do |ext|
      prefix = prefix.split(ext)[0]
    end
    ext = File.extname(filename)
    prefix = prefix.split(ext)[0] unless ext.empty?
    prefix
  end

  def link_original(filename)
    original = nil
    if File.exist?(filename)
      links = `ls -l #{filename}`.chomp.split('->')
      original = links[1].strip if links.size > 1
    end
    original
  end

  def l2name(l)
    lname = %w(s p d f g h i j)
    lname[l]
  end
end
