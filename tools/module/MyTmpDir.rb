#!/usr/bin/env ruby

#
# my tmpdir on bern series
#   written by Hiroki Funashima, Kobe University, 2016
#
module MyTmpDir
  require 'fileutils'
  include  FileUtils
  def my_tmp_dir
    "/Users/Shared/#{ENV['USER']}" + ENV['PWD'].split(ENV['HOME'])[1]
  end

  def generate_my_tmp_dir
    mkdir_p(my_tmp_dir)
  end
end
