#!/usr/bin/env ruby
#
#
module KsMailUtil
  require 'mail'
  def send_result(configure, from_add, to_add, logfile = nil)
    pwd = ENV['PWD'].split(ENV['HOME'])[1]
    series = configure[:prefix] + configure[:index].to_s
    logfile = series + '.log' if logfile.nil?
    mail = Mail.new
    mail.from = from_add
    mail.to = to_add
    mail.subject = "end of calculation: #{series} at ~#{pwd}"
    data = File.read(logfile) if File.exist?(logfile)
    dstfile = configure[:dir] + '/' + series + 'dst.dat'
    data << File.read(dstfile) if File.exist?(dstfile)
    mail.body = data
    mail.deliver
  end
end
