#!/usr/bin/env ruby2.3
#
#
module ResultMail
    def send_result(configure)
        to_add = ENV['USER']
        result_file = configure[:dir] + '/' + configure[:prefix] + configure[:index] + 'dst.dat'
        subject = "'End of calculation at #{ENV['PWD']}'"
        `cat #{result_file} | mail -s #{subject} #{to_add}`
    end
end
